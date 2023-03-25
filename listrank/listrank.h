#include "parallel.h"
#include "random.h"

// For timing parts of your code.
#include "get_time.h"

// For computing sqrt(n)
#include <math.h>

#include <thread>

using namespace parlay;

// Some helpful utilities
namespace {

// returns the log base 2 rounded up (works on ints or longs or unsigned
// versions)
template <class T>
size_t log2_up(T i) {
  assert(i > 0);
  size_t a = 0;
  T b = i - 1;
  while (b > 0) {
    b = b >> 1;
    a++;
  }
  return a;
}

}  // namespace

struct ListNode {
  ListNode* next;
  size_t rank;
  ListNode(ListNode* next) : next(next), rank(std::numeric_limits<size_t>::max()) {}
};

// Serial List Ranking. The rank of a node is its distance from the
// tail of the list. The tail is the node with `next` field nullptr.
//
// The work/depth bounds are:
// Work = O(n)
// Depth = O(n)
void SerialListRanking(ListNode* head) {
  size_t ctr = 0;
  ListNode* save = head;
  while (head != nullptr) {
    head = head->next;
    ++ctr;
  }
  head = save;
  --ctr;  // last node is distance 0
  while (head != nullptr) {
    head->rank = ctr;
    head = head->next;
    --ctr;
  }
}

// Wyllie's List Ranking. Based on pointer-jumping.
//
// The work/depth bounds of your implementation should be:
// Work = O(n*\log(n))
// Depth = O(\log^2(n))
void WyllieListRanking(ListNode* L, size_t n) {
  size_t* succ = (size_t*)malloc(n * sizeof(size_t));
  size_t* dist = (size_t*)malloc(n * sizeof(size_t));

  // populate the succ array
  auto f1 = [&](size_t i) {
    if (L[i].next == nullptr) {
      succ[i] = n;
    } else {
      succ[i] = L[i].next - L;
    }
  };

  parallel_for(0, n, f1);

  // populate the dist array
  auto f2 = [&](size_t i) {
    if (L[i].next == nullptr) {
      dist[i] = 0;
    } else {
      dist[i] = 1;
    }
  };

  parallel_for(0, n, f2);

  size_t* succ_prime = (size_t*)malloc(n * sizeof(size_t));
  size_t* dist_prime = (size_t*)malloc(n * sizeof(size_t));

  auto f3 = [&](size_t i) {
    if (succ[i] == n) {
      dist_prime[i] = dist[i];
      succ_prime[i] = n;
    } else {
      dist_prime[i] = dist[i] + dist[succ[i]];
      succ_prime[i] = succ[succ[i]];
    }
  };

  for (size_t k=0; k < log2_up(n); k++) {
    parallel_for(0, n, f3);
    std::swap(succ, succ_prime);
    std::swap(dist, dist_prime);
  }

  auto f4 = [&](size_t i) {
    L[i].rank = dist[i];
  };

  parallel_for(0, n, f4);

  free(succ);
  free(succ_prime);
  free(dist);
  free(dist_prime);
}


struct LinkedListNode {
  LinkedListNode* next;
  size_t rank;
  size_t weight;
  bool flag;
};

// Sampling-Based List Ranking
//
// The work/depth bounds of your implementation should be:
// Work = O(n) whp
// Depth = O(\sqrt(n)* \log(n)) whp
void SamplingBasedListRanking(ListNode* L, size_t n, long num_samples=-1, parlay::random r=parlay::random(0)) {
  // auto num_threads = std::thread::hardware_concurrency();
  // std::cout << num_threads << std::endl;
  
  // Perhaps use a serial base case for small enough inputs?

  if (num_samples == -1) {
    num_samples = sqrt(n);
  }

  // Generate 0 with probability (n-sqrt(n)) / n and 1 with probability (sqrt(n)/n)

  // std::random_device rd;
  // std::default_random_engine e1(rd());
  // std::discrete_distribution<size_t> discrete_dist({n - double(num_samples), double(num_samples)});

  // Create an array to store each LinkedListNode
  LinkedListNode* sampled_list = (LinkedListNode*)malloc(n * sizeof(LinkedListNode));
  size_t head_pos = -1;
  
  // true flag means that node is sampled. false flag means otherwise.
  auto f1 = [&](size_t i) {
    if (L[i].next == nullptr) { // we've found the tail of the list. add it to our sampled nodes
      sampled_list[i].flag = true;
    } else if (L == &L[i]) { // we've found the head of our list
      head_pos = i;
      sampled_list[i].flag = true;
    } else {
      size_t rand_val = r.ith_rand(i) % num_samples;

      if (rand_val == 0) {
        sampled_list[i].flag = true;
      } else {
        sampled_list[i].flag = false;
      }
    }
  };

  parallel_for(0, n, f1);

  // Assign the weights and next pointers for each of the sampled nodes

  auto assign_weights = [&](size_t i) {
    if (sampled_list[i].flag) {
      ListNode* sample_head = &L[i];
      size_t num_nonsampled = 0;

      while (sample_head->next != nullptr && !sampled_list[sample_head->next - L].flag) {
        num_nonsampled += 1;
        sample_head = sample_head->next;
      }

      sampled_list[i].weight = num_nonsampled;
      
      if (sample_head->next == nullptr) {
        sampled_list[i].next = nullptr;
      } else {
        sampled_list[i].next = &sampled_list[sample_head->next - L];
      }
    }
  };

  parallel_for(0, n, assign_weights);

  // Run the serial list-ranking algo on the sample only list

  LinkedListNode* head = &sampled_list[head_pos];

  size_t ctr = 0;
  LinkedListNode* save = head;
  while (head != nullptr) {
    ctr = ctr + head->weight + 1;
    head = head->next;
  }
  head = save;
  --ctr;  // last node is distance 0
  while (head != nullptr) {
    head->rank = ctr;
    ctr = ctr - head->weight - 1;
    head = head->next;
  }

  // Assign the ranks for the remaining non-sampled nodes

  auto assign_ranks = [&](size_t i) {
    if (sampled_list[i].flag) {
      ListNode* sample_head = &L[i];
      size_t num_nonsampled = 0;

      sample_head->rank = sampled_list[i].rank;

      while (sample_head->next != nullptr && !sampled_list[sample_head->next - L].flag) {
        num_nonsampled += 1;
        sample_head = sample_head->next;
        sample_head->rank = sampled_list[i].rank - num_nonsampled;
      }
    }
  };

  parlay::parallel_for(0, n, assign_ranks);
}

