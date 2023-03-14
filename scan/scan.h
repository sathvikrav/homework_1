#include "parallel.h"

using namespace parlay;

// A serial implementation for checking correctness.
//
// Work = O(n)
// Depth = O(n)
template <class T, class F>
T scan_inplace_serial(T *A, size_t n, const F& f, T id) {
  T cur = id;
  for (size_t i=0; i<n; ++i) {
    T next = f(cur, A[i]);
    A[i] = cur;
    cur = next;
  }
  return cur;
}

// Parallel in-place prefix sums. Your implementation can allocate and
// use an extra n*sizeof(T) bytes of memory.
//
// The work/depth bounds of your implementation should be:
// Work = O(n)
// Depth = O(\log(n))

template <class T, class F>
T scan_inplace(T *A, size_t n, const F& f, T id) {
  // return id;  // TODO
  // return id;
  T* L = (T*)malloc((n - 1) * sizeof(T));
  T total = scan_up(A, n, L, n - 1, f);
  scan_down(A, n, L, n - 1, f, id);
  return total;
}

template <class T, class F>
size_t scan_up(T *A, size_t n, T *L, size_t k, const F& f) {
  if (n == 1) {
    return A[0];
  } else {
    size_t m = n / 2;
    T l, r;
    auto f1 = [&]() { l = scan_up(A, m, L, m - 1, f); };
    auto f2 = [&]() { r = scan_up(A + m, n, L + m, k, f); };
    par_do(f1, f2);
    L[m - 1] = l;
    return f(l, r);
  }
}

template <class T, class F>
void scan_down(T *A, size_t n, T *L, size_t k, const F& f, T id) {
  if (n == 1) {
    A[0] = id;
    return;
  } else {
    size_t m = n / 2;
    auto f1 = [&]() { scan_down(A, m, L, m - 1, f, id); };
    auto f2 = [&]() { scan_down(A + m, n, L + m, k, f, f(id, L[m - 1])); };
    par_do(f1, f2);
    return;
  }
}
