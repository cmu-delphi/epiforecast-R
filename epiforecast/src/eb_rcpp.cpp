// author_header begin
// Copyright (C) 2017 Logan C. Brooks
//
// This file is part of epiforecast.  Algorithms included in epiforecast were developed by Logan C. Brooks, David C. Farrow, Sangwon Hyun, Shannon Gallagher, Ryan J. Tibshirani, Roni Rosenfeld, and Rob Tibshirani (Stanford University), members of the Delphi group at Carnegie Mellon University.
//
// Research reported in this publication was supported by the National Institute Of General Medical Sciences of the National Institutes of Health under Award Number U54 GM088491. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health. This material is based upon work supported by the National Science Foundation Graduate Research Fellowship Program under Grant No. DGE-1252522. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundation. David C. Farrow was a predoctoral trainee supported by NIH T32 training grant T32 EB009403 as part of the HHMI-NIBIB Interfaces Initiative. Ryan J. Tibshirani was supported by NSF grant DMS-1309174.
// author_header end
// license_header begin
// epiforecast is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 2 of the License.
//
// epiforecast is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with epiforecast.  If not, see <http://www.gnu.org/licenses/>.
// license_header end

#include <iterator>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <array>
#include <cassert>
#include <limits>
#include <type_traits>
#include <functional>
#include <map>
#include <cmath>

#include <Rcpp.h>

#include "epiforecast_types.hpp"

#define OUTELT(elt) \
  do { \
    std::cout << #elt << ": " << elt << std::endl;	\
  }while(0)
#define OUTITBL(itbl) \
  do { \
    std::cout << #itbl << ": "; \
    for (auto &&elt : itbl) { \
      std::cout << elt << ", "; \
    } \
    std::cout << std::endl; \
  }while(0)

// xxx remove
// auto &cout = std::cout;

// dims_diff(t0, t1) returns the difference between the two iterators, like
// t0-t1 or std::distance(t1, t0), but along each dimension; the result should
// be an array of std::ptrdiff_t. Overloads are done by providing a template
// specialization or static member function.
// template <typename Iterator>
// auto dims_diff(const Iterator &it0, const Iterator &it1) -> decltype(it0.dims_diff(it1)) {
//   return it0.dims_diff(it1);
// }
// // dims_diff for std::vector iterators:
// template <typename Elt, typename Allocator>
// auto dims_diff(const typename std::vector<Elt, Allocator>::iterator &v0, const typename std::vector<Elt, Allocator>::iterator &v1) -> decltype(v0 - v1) {
//   return v0 - v1;
// }
// // dims_diff for std::vector const iterators:
// template <typename Elt, typename Allocator>
// auto dims_diff(const typename std::vector<Elt, Allocator>::const_iterator &v0, const typename std::vector<Elt, Allocator>::const_iterator &v1) -> decltype(v0 - v1) {
//   return v0 - v1;
// }

// todo standardize dims_diff, ndims, etc.

// NDims<IterableOrRef>::value is the dimensionality of IterableOrRef; it only
// handles types with dimensionality that can be determined at compile time. It
// forwards the cv&ref-stripped typename to NDimsHelper. Overloads are done by
// providing a template specialization or static member function.
template <typename Iterable>
struct NDimsHelper {
  static constexpr std::size_t value = Iterable::ndims();
};
template <typename Elt, typename Allocator>
struct NDimsHelper<std::vector<Elt, Allocator> > {
  static constexpr std::size_t value = 1;
};
template<typename IterableOrRef>
struct NDims {
  static constexpr std::size_t value = NDimsHelper<typename std::remove_cv<typename std::remove_reference<IterableOrRef>::type>::type>::value;
};
// // shorthand for NDims<IterableOrRef>::value
// template <typename IterableOrRef>
// constexpr std::size_t t_ndims = NDims<IterableOrRef>::value;

// xxx vs static non-function members, since this doesn't change
template <typename Iterable>
constexpr std::size_t ndims(const Iterable &t) {
  return t.ndims();
}
template <typename Elt, typename Allocator>
constexpr std::size_t ndims(const std::vector<Elt, Allocator> &) {
  return 1;
}

// dims(itbl) returns the size of itbl along each of its dimensions. Overloads
// are done by providing a template specialization or member function.
template <typename Iterable>
std::array<std::size_t, NDims<Iterable>::value> dims(const Iterable &t) {
  return t.dims();
}
template <typename Elt, typename Allocator>
std::array<std::size_t, 1> dims(const std::vector<Elt, Allocator> &vec) {
  return std::array<std::size_t,1>{{vec.size()}};
}

// FunctionTraits provides information about a function type (not other callable
// types): the types of its arguments and result. It probably fails or has
// strange behavior for functions with default args.
template<typename Function, typename = typename std::enable_if<std::is_function<Function>::value>::type>
struct FunctionTraits;
// xxx vs static_assert
template<typename Result, typename... Args>
struct FunctionTraits<Result(Args...)> {
  using result = Result;
  template <std::size_t arg_i>
  using arg = typename std::tuple_element<arg_i, std::tuple<Args...>>::type;
};

// `TypePrinter<ComplexTypeExpressions...> arbitrary_name;` will cause
// compilation to fail and print out what ComplexTypeExpressions are; it is
// useful for debugging.
template <typename... Types> struct TypePrinter;

// Ignore(args...) does nothing; it is typically used to either shorten code
// that throws away the results of a function call, or to make it more explicit.
template <typename... Ts> void Ignore(Ts &&...) {}

// Iterators in this file provide a proxy typedef, to which the results of
// operator* can be converted. It is basically `reference` in the standard
// iterator framework, but is not required to be an actual reference. The
// following is a member detector that produces the `proxy` type for an Iterator
// type if it is present, and otherwise gives the `reference` type.
template <typename Iterator> struct proxy_for_helper {
protected:
  template <typename IteratorAgain>
  static typename IteratorAgain::proxy ProxyOverload(typename IteratorAgain::proxy *);
  template <typename IteratorAgain>
  // static typename IteratorAgain::value_type &ProxyOverload(...); // xxx breaks on const vector<T>; gives nonconst reference
  static typename IteratorAgain::reference ProxyOverload(...);

public:
  using type = decltype(ProxyOverload<Iterator>(nullptr));
};

// shorthand
template <typename Iterator>
using proxy_for = typename proxy_for_helper<Iterator>::type;

// iterator_for<IterableOrRef> is the iterator or const_iterator type
// corresponding to the iterable (aka range), selected based on the const
// qualifier of IterableOrRef.
template <typename IterableOrRef> struct iterator_for_helper {
protected:
  using Iterable = typename std::remove_reference<IterableOrRef>::type;
  static constexpr bool is_const = std::is_const<Iterable>::value;

public:
  using type =
      typename std::conditional<is_const, typename Iterable::const_iterator,
                                typename Iterable::iterator>::type;
};
template <typename IterableOrRef>
using iterator_for = typename iterator_for_helper<IterableOrRef>::type;
template <typename IterableOrRef> struct const_iterator_for_helper {
protected:
  using Iterable = typename std::remove_reference<IterableOrRef>::type;

public:
  using type = typename Iterable::const_iterator;
};
template <typename IterableOrRef>
using const_iterator_for =
    typename const_iterator_for_helper<IterableOrRef>::type;

// todo IterablePointer iterable class template

// todo copy in more complete iterable code, stl compliance

// ShallowCacheIterator caches the results of another iterator in a member
// variable so they are not recomputed repeatedly. Since it returns a reference
// to this member, it may be helpful for making iterables that return
// non-reference types more STL-compliant. ShallowCacheIterable produces
// ShallowCacheIterator's for another Iterable. Note that it uses value
// semantics, so large iterables should be wrapped in some sort of iterable
// pointer class or moved into the ShallowCacheIterable.
// MakeShallowCacheIterable is an object generator pattern (shorthand to avoid
// explicitly stating template args).
template <typename Iterator> class ShallowCacheIterator {
  static_assert(!std::is_reference<proxy_for<Iterator>>::value,
                "ShallowCacheIterators for reference proxies are currently "
                "unimplemented.");

public:
  using value_type = typename Iterator::value_type;
  using proxy = proxy_for<Iterator> &;
  using reference = proxy &&; // (= proxy)
  // xxx make compliant, extend operations, inherit category
  using iterator_category = std::input_iterator_tag;
  using pointer = value_type *;
  using difference_type = std::ptrdiff_t;

protected:
  Iterator iterator_;
  Iterator last_deref_it_;
  typename std::remove_const<proxy_for<Iterator>>::type cache_;
  // assignable_for<proxy_for<Iterator>> cache_;
  // xxx go back to bool valid_; approach

public:
  explicit ShallowCacheIterator(const Iterator &iterator, const Iterator &end)
      : iterator_(iterator), last_deref_it_(end), cache_() {}

  proxy operator*() {
    if (iterator_ != last_deref_it_) {
      cache_ = *iterator_;
      last_deref_it_ = iterator_;
    }
    return cache_;
  }

  ShallowCacheIterator &operator++() {
    ++iterator_;
    return *this;
  }

  bool operator==(const ShallowCacheIterator &other) const {
    return iterator_ == other.iterator_;
  }
  bool operator!=(const ShallowCacheIterator &other) const {
    return !operator==(other);
  }

  difference_type operator-(const ShallowCacheIterator &other) const {
    return iterator_ - other.iterator_;
  }
  // auto dims_diff(const ShallowCacheIterator &other) const {
  //   return ::dims_diff(iterator_, other.iterator_);
  // }
};

template <typename Iterable> class ShallowCacheIterable {
public:
  using iterator = ShallowCacheIterator<iterator_for<Iterable>>;
  using const_iterator = ShallowCacheIterator<const_iterator_for<Iterable>>;

protected:
  Iterable iterable_;

public:
  template <typename IterableOrRef,
            typename = typename std::enable_if<
                std::is_same<Iterable, typename std::remove_reference<
                                           IterableOrRef>::type>::value>::type>
  explicit ShallowCacheIterable(IterableOrRef &&iterable)
      : iterable_(std::forward<IterableOrRef>(iterable)) {}
  ShallowCacheIterable(const ShallowCacheIterable &other) = default;
  ShallowCacheIterable(ShallowCacheIterable &&other) = default;
  ShallowCacheIterable &operator=(const ShallowCacheIterable &other) = default;
  ShallowCacheIterable &operator=(ShallowCacheIterable &&other) = default;

  iterator begin() { return iterator(iterable_.begin(), iterable_.end()); }
  iterator end() { return iterator(iterable_.end(), iterable_.end()); }
  const_iterator cbegin() const { return const_iterator(iterable_.cbegin(), iterable_.cend()); }
  const_iterator cend() const { return const_iterator(iterable_.cend(), iterable_.cend()); }
  const_iterator begin() const { return cbegin(); }
  const_iterator end() const { return cend(); }

  std::size_t size() const { return iterable_.size(); }
  // static constexpr std::size_t ndims() { return ::ndims<Iterable>(); }
  static constexpr std::size_t ndims() { return NDims<Iterable>::value; }
  std::array<std::size_t, ndims()> dims() const {
    return ::dims(iterable_);
  }
};

template <typename IterableOrRef>
ShallowCacheIterable<typename std::remove_reference<IterableOrRef>::type>
MakeShallowCacheIterable(IterableOrRef &&iterable) {
  return ShallowCacheIterable<
      typename std::remove_reference<IterableOrRef>::type>(
      std::forward<IterableOrRef>(iterable));
}


// template <typename IteratorRange>
// class IteratorRangeCRTP {
// protected:
//   IteratorRange *nonconst_crtp_this() {
//     return static_cast<IteratorRange *>(this);
//   }
//   const IteratorRange *const_crtp_this() const {
//     return static_cast<const IteratorRange *>(this);
//   }
//   IteratorRange *crtp_this() { return nonconst_crtp_this(); }
//   const IteratorRange *crtp_this() const { return const_crtp_this(); }

// public:
//   using iterator = typename iterator_for<IteratorRange>::type;
//   std::size_t size() const {
//     return crtp_this().end() - crtp_this().begin();
//   }
// };

template <size_t it_ind> struct PreincrementHelperHelper {
  template <typename IterableTuple, typename IteratorTuple>
  PreincrementHelperHelper(IterableTuple &iterables, IteratorTuple &iterators) {
    // cout << "incrementing ind " << it_ind << " itor\n";
    auto &current_it = ++std::get<it_ind>(iterators);
    if (current_it == std::get<it_ind>(iterables).end()) {
      current_it = std::get<it_ind>(iterables).begin();
      PreincrementHelperHelper<it_ind - 1>(iterables, iterators);
    }
  }
};
template <> struct PreincrementHelperHelper<0> {
  template <typename IterableTuple, typename IteratorTuple>
  PreincrementHelperHelper(IterableTuple &, IteratorTuple &iterators) {
    // cout << "incrementing ind 0 itor\n";
    ++std::get<0>(iterators);
  }
};

// emulate/duplicate some C++14 features:

// enable_if_t<...>: shorthand for std::enable_if<...>::type;
template <bool enable, typename TrueType = void>
using enable_if_t = typename std::enable_if<enable, TrueType>::type;

template <typename Type>
using result_of_t = typename std::result_of<Type>::type;

// INVOKE: C++11 doesn't offer INVOKE, but essentially uses INVOKE to define
// result_of, which it does offer. Even though this is a bit circular, use
// result_of to check the outputted types of INVOKE.

// The common case: call a vanilla function/ref on vanilla args/refs.
template <typename FunctionOrRef, typename... ArgsOrRefs>
inline result_of_t<FunctionOrRef&&(ArgsOrRefs&&...)>
INVOKE(FunctionOrRef &&f, ArgsOrRefs &&... args_or_refs) {
  return std::forward<FunctionOrRef>(f)(std::forward<ArgsOrRefs>(args_or_refs)...);
}

// pointer to member function:
template <typename PointerToMemberFunction, class Base, typename DerivedOrRef, typename... ArgsOrRefs,
          typename = enable_if_t<std::is_base_of<Base, std::decay_t<DerivedOrRef>>::value>>
inline auto
INVOKE(PointerToMemberFunction Base::*member_function_ptr, DerivedOrRef &&lhs_object_or_ref, ArgsOrRefs &&... rhs_args_or_refs)
  -> decltype((std::forward<DerivedOrRef>(lhs_object_or_ref).*member_function_ptr)(std::forward<ArgsOrRefs>(rhs_args_or_refs)...)) {
  return (std::forward<DerivedOrRef>(lhs_object_or_ref).*member_function_ptr)(std::forward<ArgsOrRefs>(rhs_args_or_refs)...);
}

// other pointer to member function cases: todo

// pointer to data member: todo

// VariadicProduct(args...) is the product of `args`, a variadic version of
// arg0*arg1.
std::size_t VariadicProduct() {
  return 1;
}
template <typename T, typename... Ts>
std::size_t VariadicProduct(const T &head, const Ts &... tail) {
  return head * VariadicProduct(tail...);
}
// xxx ... T & Ts std::size_t...

template<typename... Iterable>
struct DifferenceHelperHelper;

template<typename Iterable>
struct DifferenceHelperHelper<Iterable> {
  static std::ptrdiff_t Calc(std::ptrdiff_t accu,
          const Iterable &itbl, const iterator_for<Iterable> &it0, const iterator_for<Iterable> &it1) {
    // cout << "lastit: " << accu << " * " << itbl.size() << " + " << (it0 - it1) << " = " << accu * static_cast<std::ptrdiff_t>(itbl.size()) + (it0 - it1) << "\n";
    return accu * static_cast<std::ptrdiff_t>(itbl.size()) + (it0 - it1);
  }
};

template<typename HeadIterable, typename... TailIterables>
struct DifferenceHelperHelper<HeadIterable, TailIterables...> {
  static std::ptrdiff_t Calc(std::ptrdiff_t accu,
     const HeadIterable &head_itbl, const TailIterables &... tail_itbls,
     const iterator_for<HeadIterable> &head_it0, const iterator_for<TailIterables> &... tail_its0,
     const iterator_for<HeadIterable> &head_it1, const iterator_for<TailIterables> &... tail_its1
          ) {
    // cout << "accu': " << accu << " * " << head_itbl.size() << " + " << (head_it0 - head_it1) << " = " << accu * static_cast<std::ptrdiff_t>(head_itbl.size()) + (head_it0 - head_it1) << "\n";
    return DifferenceHelperHelper<TailIterables...>::Calc(accu * static_cast<std::ptrdiff_t>(head_itbl.size()) + (head_it0 - head_it1), tail_itbls..., tail_its0..., tail_its1...);
  }
};

// MakeIter(callable).iter(args) calls `callable` on each of `args`, intended
// for a `callable` that has a void return type. It is analogous to `List.iter`
// and `Array.iter` in ML.
template <typename Callable>
class Iter {
protected:
  Callable callable_;
public:
  Iter(const Callable &callable): callable_(callable) {}

  void iter() const {
  }

  template<typename Head, typename... Tail>
  void iter(Head &&head, Tail &&... tail) {
    callable_(std::forward<Head>(head));
    iter(tail...);
  }
};

template <typename Callable>
Iter<Callable> MakeIter(const Callable &callable) {
  return Iter<Callable>(callable);
}

// xxx not all callable types currently supported; e.g., functions must be made
// into function references or pointers

// MakeCartesianProductIterable(callable, iterables...) maps `callable` over the
// cartesian product of `iterables`. If `callable` is `std::make_tuple`, then
// this is just the cartesian product of `iterables`.
template <bool is_const, typename Callable, typename... Iterables>
class CartesianProductIterator {
// protected:
public:
  // xxx vs forward-declaring CartesianProductIterable, using its types, or taking in CartesianProductIterable as a template arg?
  // const Callable &callable_;
  // Callable &callable_;
  // Callable &callable_;
  // xxx vs std::function
  Callable callable_;
  typename std::conditional<is_const, const std::tuple<Iterables...>, std::tuple<Iterables...>>::type &iterables_;
  // using temp = std::tuple<iterator_for<typename std::conditional<is_const, const Iterables, Iterables>::type>...>;
  std::tuple<typename std::conditional<is_const, const_iterator_for<Iterables>, iterator_for<Iterables> >::type...> iterators_;

public:
  // using value_type = std::tuple<typename
  // iterator_for<Iterables>::value_type...>;
  // using proxy = std::tuple<proxy_for<iterator_for<Iterables>>...>;
  using proxy = typename std::result_of<Callable(
      proxy_for<iterator_for<Iterables>>...)>::type;
  using value_type = proxy; // xxx revisit; if callable specifies result type
                            // (e.g., functor typedef), use that instead
  using reference = proxy &&;
  // xxx make compliant, extend operations, inherit category
  using iterator_category = std::input_iterator_tag;
  using pointer = value_type*;
  using difference_type = std::ptrdiff_t;

protected:
  template <size_t... it_inds>
  proxy DereferenceHelper(std::index_sequence<it_inds...>) {
    return INVOKE(callable_, std::forward<proxy_for<iterator_for<Iterables>>>(
                                 *std::get<it_inds>(iterators_))...);
  }
  template <size_t it_ind> void PreincrementHelper() {
    PreincrementHelperHelper<it_ind>(iterables_, iterators_);
  }
  template <size_t... it_inds>
  std::ptrdiff_t DifferenceHelper(const CartesianProductIterator &other, std::index_sequence<it_inds...>) const {
    return DifferenceHelperHelper<Iterables...>::Calc(0, std::get<it_inds>(iterables_)..., std::get<it_inds>(iterators_)..., std::get<it_inds>(other.iterators_)...);
  }
  // template <size_t... it_inds>
  // auto
  // dims_diff_helper(const CartesianProductIterator &other, std::index_sequence<it_inds...>) const {
  //   return std::vector<std::ptrdiff_t>({{std::get<it_inds>(iterators_)-std::get<it_inds>(other.iterators_)...}});
  // }
  // xxx default out the index_sequences on the helpers

public:
  explicit CartesianProductIterator(
            // Callable &callable,
            const Callable &callable,
            decltype(iterables_) iterables,
            const decltype(iterators_) &iterators)
      : callable_(callable), iterables_(iterables), iterators_(iterators) {
    static_assert(std::is_reference<decltype(iterables)>::value, "iterables should have been a reference");
  }
  CartesianProductIterator(const CartesianProductIterator &) = default;
  CartesianProductIterator(CartesianProductIterator &&) = default;
  CartesianProductIterator &operator=(const CartesianProductIterator &other) {
    assert(&iterables_ == &other.iterables_);
    iterators_ = other.iterators_;
    return *this;
  }
  CartesianProductIterator &operator=(CartesianProductIterator &&other) {
    assert(&iterables_ == &other.iterables_);
    iterators_ = std::move(other.iterators_);
    return *this;
  }

  proxy operator*() {
    return DereferenceHelper(std::index_sequence_for<Iterables...>());
  }

  CartesianProductIterator &operator++() {
    static_assert(sizeof...(Iterables) != 0, "failed: sizeof...(Iterables) != 0");
    PreincrementHelper<sizeof...(Iterables)-1>();
    return *this;
  }

  difference_type operator-(const CartesianProductIterator &other) const {
    return DifferenceHelper(other, std::index_sequence_for<Iterables...>());
  }
  // auto dims_diff(const CartesianProductIterator &other) const {
  //   return dims_diff_helper(other, std::index_sequence_for<Iterables...>());
  // }

  bool operator==(const CartesianProductIterator &other) const {
    assert(&iterables_ == &other.iterables_);
    return iterators_ == other.iterators_;
  }
  bool operator!=(const CartesianProductIterator &other) const {
    return !operator==(other);
  }
};

// TemplateSum<sizes...> is the sum of `sizes`.
template <std::size_t...>
struct TemplateSum;

template <>
struct TemplateSum<> {
  static constexpr std::size_t value = 0;
};
template <std::size_t head, std::size_t... tail>
struct TemplateSum<head, tail...> {
  static constexpr std::size_t value = head + TemplateSum<tail...>::value;
};

template <typename Callable, typename... Iterables>
class CartesianProductIterable {
  // xxx turn into enable_if?
  static_assert(!std::is_reference<Callable>::value, "callable references not supported");
  static_assert(sizeof...(Iterables) != 0, "CartesianProductIterable requires >=1 iterable.");

protected:
  Callable callable_;
  std::tuple<Iterables...> iterables_;

public:
  using iterator = CartesianProductIterator<false, Callable, Iterables...>;
  // using const_iterator = CartesianProductIterator<true, const Callable, const Iterables...>;
  // using const_iterator = CartesianProductIterator<true, const Callable, Iterables...>;
  using const_iterator = CartesianProductIterator<true, Callable, Iterables...>;

protected:
  template <size_t... it_inds>
  std::tuple<iterator_for<Iterables>...>
      BeginHelper(std::index_sequence<it_inds...>) {
    return std::tuple<iterator_for<Iterables>...>(
        std::get<it_inds>(iterables_).begin()...);
  }
  template <size_t... it_inds>
  std::tuple<const_iterator_for<Iterables>...>
      CBeginHelper(std::index_sequence<it_inds...>) const {
    return std::tuple<const_iterator_for<Iterables>...>(
        std::get<it_inds>(iterables_).cbegin()...);
  }

  template <size_t... it_inds>
  std::size_t SizeHelper(std::index_sequence<it_inds...>) const {
    return VariadicProduct(std::get<it_inds>(iterables_).size()...);
  }

public:
  // xxx checks that same iterables
  template <typename... IterablesOrRefs>
  explicit CartesianProductIterable(const Callable &callable,
                                    IterablesOrRefs &&... iterables)
      : callable_(callable),
        iterables_(std::forward<IterablesOrRefs>(iterables)...) {
  }
  CartesianProductIterable(const CartesianProductIterable &other) = default;
  CartesianProductIterable(CartesianProductIterable &&other) = default;
  CartesianProductIterable &operator=(const CartesianProductIterable &other)
  = default;
  CartesianProductIterable &operator=(CartesianProductIterable &&other) =
  default;

  iterator begin() {
    return iterator(callable_, iterables_,
                    BeginHelper(std::index_sequence_for<Iterables...>()));
  }
  iterator end() {
    auto its = BeginHelper(std::index_sequence_for<Iterables...>());
    if (sizeof...(Iterables) != 0) {
      std::get<0>(its) = std::get<0>(iterables_).end();
    }
    return iterator(callable_, iterables_, its);
  }
  const_iterator cbegin() const {
    const_iterator result(callable_, iterables_, CBeginHelper(std::index_sequence_for<Iterables...>()));
    return result;
  }
  const_iterator cend() const {
    auto its = CBeginHelper(std::index_sequence_for<Iterables...>());
    std::get<0>(its) = std::get<0>(iterables_).cend();
    return const_iterator(callable_, iterables_, its);
  }
  const_iterator begin() const {
    return cbegin();
  }
  const_iterator end() const {
    return cend();
  }

  std::size_t size() const {
    return SizeHelper(std::index_sequence_for<Iterables...>());
  }
  static constexpr std::size_t ndims() { return TemplateSum<NDims<Iterables>::value...>::value; }
protected:
  template<std::size_t... it_inds>
  std::array<std::size_t, ndims()> DimsHelper(std::index_sequence<it_inds...>) const {
    std::array<std::size_t, ndims()> result;
    auto result_it = result.begin();
    auto appender_lambda = [&result_it](const auto &itbl_dims) -> void {
      // xxx std::copy
  for (const std::size_t &dim : itbl_dims) {
    *result_it++ = dim;
  }
    };
    // appender_lambda(std::array<std::size_t,1>{{1}});
    // appender_lambda(std::array<std::size_t,2>{{5,2}});
    MakeIter(appender_lambda).iter(::dims(std::get<it_inds>(iterables_))...);
    return result;
  }
public:
  std::array<std::size_t, ndims()> dims() const {
    return DimsHelper(std::index_sequence_for<Iterables...>());
  }
};

template <typename... IterablesOrRefs>
CartesianProductIterable<
    decltype(
        &std::forward_as_tuple<proxy_for<iterator_for<IterablesOrRefs>>...>),
    typename std::remove_reference<IterablesOrRefs>::type...>
MakeCartesianProductIterable(IterablesOrRefs &&... iterables) {
  return CartesianProductIterable<
      decltype(
          &std::forward_as_tuple<proxy_for<iterator_for<IterablesOrRefs>>...>),
      typename std::remove_reference<IterablesOrRefs>::type...>(
      &std::forward_as_tuple<proxy_for<iterator_for<IterablesOrRefs>>...>,
      std::forward<IterablesOrRefs>(iterables)...);
}

template <typename Callable, typename... IterablesOrRefs>
CartesianProductIterable<
    Callable, typename std::remove_reference<IterablesOrRefs>::type...>
MakeFMapCartesianProductIterable(const Callable &callable,
                                 IterablesOrRefs &&... iterables) {
  return CartesianProductIterable<
      Callable, typename std::remove_reference<IterablesOrRefs>::type...>(
      callable, std::forward<IterablesOrRefs>(iterables)...);
}

// Iota<>(n) is an iterable over the first n nonnegative integers, similar to
// `seq_len` in R but does not actually construct an n-length vector.
template <typename Output=std::size_t> class IotaIterator {
public:
  using value_type = const Output;
  using proxy = const Output;
  using reference = proxy &&;
  // xxx make compliant, extend operations, inherit category
  using iterator_category = std::input_iterator_tag;
  using pointer = value_type*;
  using difference_type = std::ptrdiff_t;

protected:
  std::size_t index_;

public:
  explicit IotaIterator(std::size_t index) : index_(index) {}

  proxy operator*() { return Output(index_); }
  IotaIterator &operator++() {
    ++index_;
    return *this;
  }
  bool operator==(const IotaIterator &other) const {
    return index_ == other.index_;
  }
  bool operator!=(const IotaIterator &other) const {
    return !operator==(other);
  }

  difference_type operator-(const IotaIterator &other) const {
    return static_cast<difference_type>(index_-other.index_);
  }
  // auto dims_diff(const IotaIterator &other) const {
  //   return operator-(other);
  // }
};

template <typename Output = std::size_t> class Iota {
public:
  using iterator = IotaIterator<Output>;
  using const_iterator = IotaIterator<Output>;

  const std::size_t end_;

  explicit Iota(const std::size_t &end) : end_(end) {}

  iterator begin() { return iterator(0); }
  iterator end() { return iterator(end_); }
  const_iterator cbegin() const { return const_iterator(0); }
  const_iterator cend() const { return const_iterator(end_); }
  const_iterator begin() const { return cbegin(); }
  const_iterator end() const { return cend(); }

  std::size_t size() const { return static_cast<std::size_t>(end_); }
  static constexpr std::size_t ndims() { return 1; }
  std::array<std::size_t, ndims()> dims() const {
    return std::array<std::size_t, ndims()>{{end_}};
  }
};

// todo sampling methods given an IO&RNG object.

// GaussianObservable's provide arithmetic over independent Gaussian RV's.
class GaussianObservable {
 // protected:
 public:
  Mean mean_;
  SD sd_;

 public:
  explicit GaussianObservable(Mean mean, SD sd) : mean_(mean), sd_(sd) {
    // xxx vs <= sd_epsilon; <= sd_epsilon breaks on math manipulations multiplying by 0, etc....; but fp stuff may happen and make it slightly less than zero (catastrophic cancellation, etc....); consider public vs. protected constructors for checking on creation only, proxied calculations to make sure no final results use?
    // if (sd <= 0) {
    //   throw std::invalid_argument("sd must be strictly positive and > sd_epsilon");
    // }
  }
  // GaussianObservable() : GaussianObservable(0,0) {}

  // xxx should change for case consistency
  const Mean &mean() { return mean_; }
  const SD &sd() { return sd_; }
  LogLikelihood log_likelihood(const Mean &observation) const {
    auto z = (observation - mean_) / sd_;
    LogLikelihood logphiz = (-z*z - log(2*M_PI))/2;
    return logphiz - log(sd_);
  }

  friend GaussianObservable operator+(const GaussianObservable &obs0,
                                      const GaussianObservable &obs1) {
    return GaussianObservable(obs0.mean_ + obs1.mean_, obs0.sd_ + obs1.sd_);
  }
  GaussianObservable &operator+=(const GaussianObservable &other) {
    mean_ += other.mean_;
    sd_ += other.sd_;
    return *this;
  }
  friend GaussianObservable operator-(const GaussianObservable &obs0,
                                      const GaussianObservable &obs1) {
    return GaussianObservable(obs0.mean_ - obs1.mean_, obs0.sd_ + obs1.sd_);
  }
  GaussianObservable &operator-=(const GaussianObservable &other) {
    mean_ -= other.mean_;
    sd_ += other.sd_;
    return *this;
  }
  friend GaussianObservable operator+(const GaussianObservable &obs,
                                      const Mean &mean) {
    return GaussianObservable(obs.mean_ + mean, obs.sd_);
  }
  friend GaussianObservable operator+(const Mean &mean,
                                      const GaussianObservable &obs) {
    return GaussianObservable(mean + obs.mean_, obs.sd_);
  }
  GaussianObservable &operator+=(const Mean &mean) {
    mean_ += mean;
    return *this;
  }
  friend GaussianObservable operator-(const Mean &mean,
                                      const GaussianObservable &obs) {
    return GaussianObservable(mean - obs.mean_, obs.sd_);
  }
  friend GaussianObservable operator-(const GaussianObservable &obs,
                                      const Mean &mean) {
    return GaussianObservable(obs.mean_ - mean, obs.sd_);
  }
  GaussianObservable &operator-=(const Mean &mean) {
    mean_ -= mean;
    return *this;
  }
  friend GaussianObservable operator*(const Mean &scale,
                                      const GaussianObservable &obs) {
    return GaussianObservable(scale * obs.mean_, fabs(scale) * obs.sd_);
  }
  friend GaussianObservable operator*(const GaussianObservable &obs,
                                      const Mean &scale) {
    return GaussianObservable(obs.mean_ * scale, obs.sd_ * fabs(scale));
  }
  GaussianObservable &operator*=(const Mean &scale) {
    mean_ *= scale;
    sd_ *= fabs(scale);
    return *this;
  }
  friend GaussianObservable operator/(const GaussianObservable &obs,
                                      const Mean &scale) {
    return GaussianObservable(obs.mean_ / scale, obs.sd_ / fabs(scale));
  }
  GaussianObservable &operator/=(const Mean &scale) {
    mean_ /= scale;
    sd_ /= fabs(scale);
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &out,
                                  const GaussianObservable &obs) {
    return out << "N<m" << obs.mean_ << ",s" << obs.sd_ << ">";
  }
  // xxx meaning of rel ops...
  friend bool operator<(const GaussianObservable &obs, const Mean &mean) {
    return obs.mean_ < mean;
  }
  friend bool operator>(const GaussianObservable &obs, const Mean &mean) {
    return obs.mean_ > mean;
  }
  friend bool operator<=(const GaussianObservable &obs, const Mean &mean) {
    return obs.mean_ <= mean;
  }
  friend bool operator>=(const GaussianObservable &obs, const Mean &mean) {
    return obs.mean_ >= mean;
  }
  friend bool operator<(const Mean &mean, const GaussianObservable &obs) {
    return mean < obs.mean_;
  }
  friend bool operator>(const Mean &mean, const GaussianObservable &obs) {
    return mean > obs.mean_;
  }
  friend bool operator<=(const Mean &mean, const GaussianObservable &obs) {
    return mean <= obs.mean_;
  }
  friend bool operator>=(const Mean &mean, const GaussianObservable &obs) {
    return mean >= obs.mean_;
  }
  // xxx meaning of rel ops...
  bool operator<(const GaussianObservable &other) const {
    return mean_ < other.mean_;
  }
  bool operator>(const GaussianObservable &other) const {
    return mean_ > other.mean_;
  }
  bool operator<=(const GaussianObservable &other) const {
    return mean_ <= other.mean_;
  }
  bool operator>=(const GaussianObservable &other) const {
    return mean_ >= other.mean_;
  }

  struct MeanComparator {
    bool operator()(const GaussianObservable& observable0, const GaussianObservable& observable1) const {
      return observable0.mean_ < observable1.mean_;
    }
  };
};

using Observable = GaussianObservable;
using Trajectory = std::vector<std::pair<Time, Observation>>;
using MeanCurve = std::vector<std::pair<Time, Mean>>;
using ObservableCurve = std::vector<std::pair<Time, Observable>>;

// PeakIt(curve) returns an iterator to an (x,y) pair in `curve` such that y is
// >= any other y' in `curve`.
iterator_for<const ObservableCurve>
PeakIt(const ObservableCurve &curve) {
  return std::max_element(curve.cbegin(), curve.cend(),
                          [](const auto &mapping0, const auto &mapping1) {
                            return mapping0.second < mapping1.second;
                          });
}

// XTransformer provides lazy x-shifting and x-scaling of curves.
class XTransformer {
protected:
  const ObservableCurve *curve_;
  Time curve_peak_time_;
  Time x_shift_;
  Rp x_scale_;

  explicit XTransformer(const ObservableCurve *curve, Time curve_peak_time,
                        Time x_shift, Rp x_scale)
      : curve_(curve), curve_peak_time_(curve_peak_time), x_shift_(x_shift),
        x_scale_(x_scale) {}
  explicit XTransformer(const ObservableCurve &curve, Time curve_peak_time,
                        Time x_shift, Rp x_scale)
      : curve_(&curve), curve_peak_time_(curve_peak_time), x_shift_(x_shift),
        x_scale_(x_scale) {}
  // todo clean up

public:
  // xxx so can put in containers and more simply cache... but requires possibly
  // invalid pointer...
  explicit XTransformer() = default;

  explicit XTransformer(const ObservableCurve &curve)
      : curve_(&curve), curve_peak_time_(PeakIt(curve)->first), x_shift_(0),
        x_scale_(1) {}

  XTransformer XShiftedForPeakTime(Time peak_time) const {
    // cout << "XShiftedForPeakTime " << *this << "; " << peak_time << "\n";
    assert(curve_ != nullptr);
    Time result_x_shift =
      std::isnan(peak_time) ? x_shift_ : (peak_time - curve_peak_time_ * x_scale_);
    return XTransformer(*curve_, curve_peak_time_, result_x_shift, x_scale_);
  }

  XTransformer XShiftedRight(Time x_shift) const {
    assert(curve_ != nullptr);
    return XTransformer(*curve_, curve_peak_time_, this->x_shift_ + x_shift,
                        x_scale_);
  }

  XTransformer XScaledWidenAroundPeakTime(Rp x_scale) const {
    assert(curve_ != nullptr);
    Time result_x_shift =
        // curve_peak_time_ - x_scale * (curve_peak_time_ - x_shift_);
      curve_peak_time_*this->x_scale_*(1 - x_scale) + x_shift_;
    Rp result_x_scale = this->x_scale_ * x_scale;
    return XTransformer(*curve_, curve_peak_time_, result_x_shift,
                        result_x_scale);
  }

  ObservableCurve TransformedCurve() const {
    ObservableCurve result;
    result.reserve(curve_->size());
    for (const auto &mapping : *curve_) {
      result.insert(
          result.end(),
          std::make_pair(mapping.first * x_scale_ + x_shift_, mapping.second));
    }
    return result;
  }

  std::ostream &Print(std::ostream &out) const {
    assert(curve_ != nullptr);
    const_iterator_for<ObservableCurve> it = curve_->cbegin(), end = curve_->cend();
    if (it == end) {
      out << "<empty>";
    } else {
      out << "[";
      out << ((it->first) * x_scale_ + x_shift_) << "=>" << (it->second);
      ++it;
      for (; it != end; ++it) {
        out << ", ";
        out << ((it->first) * x_scale_ + x_shift_) << "=>" << (it->second);
      }
      out << "]";
    }
    return out;
  }

  friend std::ostream &operator<<(std::ostream &out,
                                  const XTransformer &x_transformer) {
    return x_transformer.Print(out);
  }
};

// cannot form pointer to constructor; use this instead for FP:
XTransformer MakeXTransformer(const ObservableCurve &curve) {
  return XTransformer(curve);
}

// xxx in zip iterator, only check first for equality, etc.

// variations on the below were tested to try to product seg faults from stack
// smashing from reference members to temporaries
auto asdf_qwer_outer_product(std::vector<int> &asdf,
                             std::vector<double> &qwer) {
  auto lambda = [](int &a, double &b) { return a * b; };
  auto outer_product_uncached =
      MakeFMapCartesianProductIterable(lambda, asdf, qwer);
  auto outer_product = MakeShallowCacheIterable(outer_product_uncached);
  return outer_product;
}

// int maino() {
//   std::vector<int> asdf = {{4, 2, 3, 1}};
//   std::vector<double> qwer = {{1.5, 2.5, 3.5}};
//   auto cart_product = MakeCartesianProductIterable(asdf, qwer);
//   for (const auto &zxcv : cart_product) {
//     std::cout << std::get<0>(zxcv) << ", " << std::get<1>(zxcv) << "\n";
//   }
//   // auto outer_product = MakeFMapCartesianProductIterable([](int &a, double
//   // &b){return a*b;},
//   //              asdf, qwer);
//   // auto outer_product =
//   // MakeShallowCacheIterable(MakeFMapCartesianProductIterable([](int &a, double
//   // &b){return a*b;}, asdf, qwer));
//   auto outer_product =
//       MakeShallowCacheIterable(MakeFMapCartesianProductIterable(
//           [](int &a, double &b) { return a * b; }, asdf, qwer));
//   for (const auto &zxcv : outer_product) {
//     std::cout << zxcv << "\n";
//   }
//   Iota<> iota(5);
//   // for (auto &zxcv : iota) {
//   // for (const auto &zxcv : iota) {
//   for (auto &&zxcv : iota) {
//     std::cout << zxcv << std::endl;
//   }
//   ShallowCacheIterable<Iota<>> caching_test(iota);
//   for (auto &zxcv : caching_test) {
//     std::cout << zxcv << std::endl;
//   }
//   std::cout << std::accumulate(outer_product.begin(), outer_product.end(), 0)
//             << std::endl;
//   auto outer_product2 = asdf_qwer_outer_product(asdf, qwer);
//   std::cout << std::accumulate(outer_product2.begin(), outer_product2.end(), 0)
//             << std::endl;
//   std::cout << "EXIT MAINO" << std::endl;
//   return 0;
// };

// CommonSize(iterables...) returns the size of each iterable in `iterables` if
// they share a common size, and throws an exception otherwise.
template<typename HeadIterable, typename... TailIterables>
struct CommonSizeHelper;

template<typename HeadIterable, typename... TailIterables>
struct CommonSizeHelper {
  std::size_t value_;
  explicit CommonSizeHelper(const HeadIterable &head, const TailIterables &... tail): value_(head.size()) {
    if (head.size() != CommonSizeHelper<TailIterables...>(tail...).value_) {
      throw std::invalid_argument("iterables do not share a common size");
    }
  }
};
template<typename Iterable>
struct CommonSizeHelper<Iterable> {
  std::size_t value_;
  explicit CommonSizeHelper(const Iterable &iterable): value_(iterable.size()) {
  }
};

// fixme check for implicit tuple const ref to temp conversion bug in ZipProductIterator instantiation

// xxx not all callable types currently supported; e.g., functions must be made
// into function references or pointers

// MakeZipProductIterable(callable, iterables...) maps callable over the "zip"
// of the iterables. If `callable` is `std::make_tuple`, then it is simply the
// zip of the iterables: an iterable whose nth element is the tuple of the nth
// element of each iterable in `iterables`.
template <typename Callable, typename... Iterables> class ZipProductIterator {
protected:
  Callable callable_;
  std::tuple<iterator_for<Iterables>...> iterators_;

public:
  // using value_type = std::tuple<typename
  // iterator_for<Iterables>::value_type...>;
  // using proxy = std::tuple<proxy_for<iterator_for<Iterables>>...>;
  using proxy = typename std::result_of<Callable(proxy_for<iterator_for<Iterables>>...)>::type;
  using value_type = proxy; // xxx revisit; if callable specifies result type
                            // (e.g., functor typedef), use that instead
  using reference = proxy &&;
  // xxx make compliant, extend operations, inherit category
  using iterator_category = std::input_iterator_tag;
  using pointer = value_type*;
  using difference_type = std::ptrdiff_t;

protected:
  template <size_t... it_inds>
  proxy DereferenceHelper(std::index_sequence<it_inds...>) {
    return INVOKE(callable_, std::forward<proxy_for<iterator_for<Iterables>>>(
                                 *std::get<it_inds>(iterators_))...);
  }
  template <size_t... it_inds>
  void PreincrementHelper(std::index_sequence<it_inds...>) {
    Ignore(++std::get<it_inds>(iterators_)...);
  }

public:
  explicit ZipProductIterator(
      const Callable &callable,
      const std::tuple<iterator_for<Iterables>...> &iterators)
      : callable_(callable), iterators_(iterators) {}
  ZipProductIterator(const ZipProductIterator &) = default;
  ZipProductIterator(ZipProductIterator &&) = default;
  ZipProductIterator &operator=(const ZipProductIterator &other) {
    // assert(&callable_ == &other.callable_);
    // assert(callable_ == other.callable_);
    iterators_ = other.iterators_;
    return *this;
  }
  ZipProductIterator &operator=(ZipProductIterator &&other) {
    // assert(&callable_ == &other.callable_);
    // assert(callable_ == other.callable_);
    iterators_ = std::move(other.iterators_);
    return *this;
  }

  proxy operator*() {
    return DereferenceHelper(std::index_sequence_for<Iterables...>());
  }

  ZipProductIterator &operator++() {
    PreincrementHelper(std::index_sequence_for<Iterables...>());
    return *this;
  }

  bool operator==(const ZipProductIterator &other) const {
    return iterators_ == other.iterators_;
  }
  bool operator!=(const ZipProductIterator &other) const {
    return !operator==(other);
  }

  difference_type operator-(const ZipProductIterator &other) const {
    return std::get<0>(iterators_) - std::get<0>(other.iterators_);
  }
  // auto dims_diff(const ZipProductIterator &other) const {
  //   return ::dims_diff(std::get<0>(iterators_), std::get<0>(other.iterators_));
  // }
};

template<typename HeadIterable, typename... TailIterables>
std::size_t CommonSize(const HeadIterable &head, const TailIterables &... tail) {
  return CommonSizeHelper<HeadIterable, TailIterables...>(head, tail...).value_;
}

template <typename Callable, typename... Iterables> class ZipProductIterable {
  // xxx turn into enable_if?
  static_assert(!std::is_reference<Callable>::value, "callable references not supported");
  static_assert(sizeof...(Iterables) != 0, "ZipProductIterable requires >=1 iterable.");

protected:
  Callable callable_;
  std::tuple<Iterables...> iterables_;

public:
  using iterator = ZipProductIterator<Callable, Iterables...>;
  using const_iterator = ZipProductIterator<Callable, const Iterables...>;

protected:
  template <size_t... it_inds>
  std::tuple<iterator_for<Iterables>...>
      BeginHelper(std::index_sequence<it_inds...>) {
    return std::tuple<iterator_for<Iterables>...>(
        std::get<it_inds>(iterables_).begin()...);
  }
  template <size_t... it_inds>
  std::tuple<iterator_for<Iterables>...>
      EndHelper(std::index_sequence<it_inds...>) {
    return std::tuple<iterator_for<Iterables>...>(
        std::get<it_inds>(iterables_).end()...);
  }
  template <size_t... it_inds>
  std::tuple<const_iterator_for<Iterables>...>
      CBeginHelper(std::index_sequence<it_inds...>) const {
    return std::tuple<const_iterator_for<Iterables>...>(
        std::get<it_inds>(iterables_).cbegin()...);
  }
  template <size_t... it_inds>
  std::tuple<const_iterator_for<Iterables>...>
      CEndHelper(std::index_sequence<it_inds...>) const {
    return std::tuple<const_iterator_for<Iterables>...>(
        std::get<it_inds>(iterables_).cend()...);
  }

public:
  // xxx checks that are same iterables
  template <typename... IterablesOrRefs>
  explicit ZipProductIterable(const Callable &callable,
                              IterablesOrRefs &&... iterables)
      : callable_(callable),
        iterables_(std::forward<IterablesOrRefs>(iterables)...) {
    Ignore(CommonSize(iterables...));
  }

  iterator begin() {
    return iterator(callable_,
                    BeginHelper(std::index_sequence_for<Iterables...>()));
  }
  iterator end() {
    return iterator(callable_,
                    EndHelper(std::index_sequence_for<Iterables...>()));
  }
  const_iterator cbegin() const {
    return const_iterator(callable_,
                    CBeginHelper(std::index_sequence_for<Iterables...>()));
  }
  const_iterator cend() const {
    return const_iterator(callable_,
                    CEndHelper(std::index_sequence_for<Iterables...>()));
  }
  const_iterator begin() const {
    return cbegin();
  }
  const_iterator end() const {
    return cend();
  }

  // xxx common size checks
  std::size_t size() const {
    return std::get<0>(iterables_).size();
  }
  static constexpr std::size_t ndims() { return NDims<typename std::tuple_element<0, decltype(iterables_)>::type>::value; }
  std::array<std::size_t, ndims()> dims() const {
    return ::dims(std::get<0>(iterables_));
  }
};

template <typename... IterablesOrRefs>
ZipProductIterable<decltype(&std::forward_as_tuple<
                            proxy_for<iterator_for<IterablesOrRefs>>...>),
                   typename std::remove_reference<IterablesOrRefs>::type...>
MakeZipProductIterable(IterablesOrRefs &&... iterables) {
  return ZipProductIterable<
      decltype(
          &std::forward_as_tuple<proxy_for<iterator_for<IterablesOrRefs>>...>),
      typename std::remove_reference<IterablesOrRefs>::type...>(
      &std::forward_as_tuple<proxy_for<iterator_for<IterablesOrRefs>>...>,
      std::forward<IterablesOrRefs>(iterables)...);
}

template <typename Callable, typename... IterablesOrRefs>
ZipProductIterable<Callable,
                   typename std::remove_reference<IterablesOrRefs>::type...>
MakeFMapZipProductIterable(const Callable &callable,
                           IterablesOrRefs &&... iterables) {
  return ZipProductIterable<
      Callable, typename std::remove_reference<IterablesOrRefs>::type...>(
      callable, std::forward<IterablesOrRefs>(iterables)...);
}

template <typename Callable, typename IterableOrRef,
    // xxx review and perhaps propagate around
    typename = typename std::result_of<Callable(proxy_for<iterator_for<IterableOrRef>>)>::type>
ZipProductIterable<Callable, typename std::remove_reference<IterableOrRef>::type>
MakeFMapIterable(const Callable &callable, IterableOrRef &&iterable) {
  return ZipProductIterable<
      Callable, typename std::remove_reference<IterableOrRef>::type>(
      callable, std::forward<IterableOrRef>(iterable));
}

// CartesianProducts and ZipProducts are interfaces to make
// CartesianProductIterables and ZipProductIterables. The `Make` methods produce
// plain cartesian products and zips; the `MakeFMap` methods map a function over
// the cartesian products and zips.
struct CartesianProducts {
  template <typename... IterablesOrRefs>
  static CartesianProductIterable<
      decltype(
          &std::forward_as_tuple<proxy_for<iterator_for<IterablesOrRefs>>...>),
      typename std::remove_reference<IterablesOrRefs>::type...>
  Make(IterablesOrRefs &&... iterables) {
    return MakeCartesianProductIterable(
        std::forward<IterablesOrRefs>(iterables)...);
  }

  template <typename Callable, typename... IterablesOrRefs>
  static CartesianProductIterable<
      Callable, typename std::remove_reference<IterablesOrRefs>::type...>
  MakeFMap(const Callable &callable, IterablesOrRefs &&... iterables) {
    return MakeFMapCartesianProductIterable(
        callable, std::forward<IterablesOrRefs>(iterables)...);
  }
};
struct ZipProducts {
  template <typename... IterablesOrRefs>
  ZipProductIterable<decltype(&std::forward_as_tuple<
                              proxy_for<iterator_for<IterablesOrRefs>>...>),
                     typename std::remove_reference<IterablesOrRefs>::
                         type...> static Make(IterablesOrRefs &&... iterables) {
    return MakeZipProductIterable(std::forward<IterablesOrRefs>(iterables)...);
  }

  template <typename Callable, typename... IterablesOrRefs>
  ZipProductIterable<
      Callable, typename std::remove_reference<IterablesOrRefs>::
                    type...> static MakeFMap(const Callable &callable,
                                             IterablesOrRefs &&... iterables) {
    return MakeFMapZipProductIterable(
        callable, std::forward<IterablesOrRefs>(iterables)...);
  }
};

// MakeObservableCurve(mean_curve, sd) produces a curve of observables, pairing
// each mean in `mean_curve` with the single standard deviation value `sd`.
ObservableCurve
MakeObservableCurve(const MeanCurve &mean_curve, SD sd) {
  auto observable_mappings = MakeFMapIterable([sd](proxy_for<iterator_for<decltype(mean_curve)>> mean_mapping) {
      return std::pair<Time,Observable>(mean_mapping.first, Observable(mean_mapping.second, sd));
  }, mean_curve);
  return ObservableCurve(observable_mappings.begin(), observable_mappings.end());
}

// ObservableCurveWithSD(curve, sd) is the result of replacing the sd's of the
// observables in `curve` with the provided `sd` unles the provided `sd` is NaN.
ObservableCurve
ObservableCurveWithSD(const ObservableCurve &curve, SD sd) {
  if (std::isnan(sd)) {
    return ObservableCurve(curve);
  } else {
    auto result_mappings = MakeFMapIterable([sd](proxy_for<iterator_for<decltype(curve)>> mapping){
    return std::make_pair(mapping.first, Observable(mapping.second.mean_, sd));
      }, curve);
    return ObservableCurve(result_mappings.begin(), result_mappings.end());
  }
}

// ObservableCurveWithSD(curve, sd) is the result of multiplying the sd's of the
// observables in `curve` with the provided `sd_scale`.
ObservableCurve
ObservableCurveWithSDScale(const ObservableCurve &curve, Rp sd_scale) {
  auto result_mappings = MakeFMapIterable([sd_scale](proxy_for<iterator_for<decltype(curve)>> mapping){
      return std::make_pair(mapping.first, Observable(mapping.second.mean_, mapping.second.sd_*sd_scale));
    }, curve);
  return ObservableCurve(result_mappings.begin(), result_mappings.end());
}

class PeakHeightScalerAboveBaseline {
protected:
  const Mean baseline_;
public:
  explicit PeakHeightScalerAboveBaseline(Mean baseline): baseline_(baseline) {}

  ObservableCurve
  operator()(const ObservableCurve &curve, Mean peak_height) const {
    if (std::isnan(peak_height)) {
      return ObservableCurve(curve);
    // } else if (peak_height < baseline_) {
    //   // xxx check on input, assert here?
    //   throw std::invalid_argument("Attempting y-scale for a peak_height below the baseline.");
    } else {
      if (peak_height < baseline_) {
        peak_height = baseline_;
      }
      Mean current_peak_height = PeakIt(curve)->second.mean_;
      Rp y_scale = (peak_height - baseline_)/(current_peak_height - baseline_);
      auto result_mappings = MakeFMapIterable([this, y_scale](proxy_for<iterator_for<decltype(curve)>> mapping){
    Observable y = mapping.second;
    Observable scaled_y = y >= baseline_? y_scale*(y-baseline_)+baseline_ : y;
    return std::make_pair(mapping.first, scaled_y);
  }, curve);
      return ObservableCurve(result_mappings.begin(), result_mappings.end());
    }
  }
};

class YScalerAboveBaseline {
protected:
  const Mean baseline_;
public:
  explicit YScalerAboveBaseline(Mean baseline): baseline_(baseline) {}

  ObservableCurve
  operator()(const ObservableCurve &curve, Rp y_scale) const {
    auto result_mappings = MakeFMapIterable([this, y_scale](proxy_for<iterator_for<decltype(curve)>> mapping){
  Observable y = mapping.second;
  Observable scaled_y = y >= baseline_? y_scale*(y-baseline_)+baseline_ : y;
  return std::make_pair(mapping.first, scaled_y);
      }, curve);
    return ObservableCurve(result_mappings.begin(), result_mappings.end());
  }
};

// MakeCeilingIterable(point_iterable, level_iterable, level_less_than_point)
// takes each point in `point_iterable` and finds the first level in
// `level_iterable` that is not less than that point, as determined by
// `level_less_than_point`.  Both iterables are assumed to be sorted.
template<typename PointIterable, typename LevelIterable, typename Callable>
auto
MakeCeilingIterable(PointIterable &point_iterable,
        LevelIterable &level_iterable,
        const Callable &level_less_than_point) {
// xxx default level_to_point: identityfn
  auto level_it = level_iterable.begin();
  // xxx http://stackoverflow.com/questions/21443023/capturing-a-reference-by-reference-in-a-c11-lambda... here and elsewhere...
  return MakeFMapIterable([&level_iterable, level_less_than_point, level_it](proxy_for<iterator_for<PointIterable>> point) mutable {
      while (level_it != level_iterable.end() && INVOKE(level_less_than_point, *level_it, point)) {
  ++level_it;
      }
      return std::make_pair(point, level_it);
    }, point_iterable);
}

// xxx checks for proper size
Observable InferWithCeiling(Time time, iterator_for<const ObservableCurve> curve_it, const ObservableCurve &curve) {
  if (curve_it == curve.begin()) {
    return curve_it->second;
  } else if (curve_it == curve.end()) {
    return (--curve_it)->second;
  } else {
    Time time1 = curve_it->first;
    const Observable &observation1 = curve_it->second;
    --curve_it;
    Time time0 = curve_it->first;
    const Observable &observation0 = curve_it->second;
    // return y0 + (x - x0) / (x1 - x0) * (y1 - y0);
    // above line causes sd expansion during interpolation
    // xxx look into more efficient way to compute the below; maybe w0 and
    // complement
    return ((time1 - time) * observation0 + (time - time0) * observation1) / (time1 - time0);
  }
}

// RepeatingIterable(elt, times) is analogous to rep(elt,times) in R, but
// without constructing a times-length vector.
template<typename Elt>
class RepeatingIterator {
public:
  using value_type = const Elt;
  using proxy = const Elt &;
  using reference = proxy;
  // xxx extend, make compliant
  using iterator_category = std::input_iterator_tag;
  using pointer = const Elt *;
  using difference_type = std::ptrdiff_t;

// xxx protected:
public:
  const Elt &elt_;
  IotaIterator<> iota_it_;

public:
  RepeatingIterator(const Elt &elt, const IotaIterator<> &iota_it): elt_(elt), iota_it_(iota_it) {
  }

  proxy operator*() {
    return elt_;
  }

  RepeatingIterator &operator++() {
    ++iota_it_;
    return *this;
  }

  bool operator==(const RepeatingIterator &other) const {
    assert(&elt_ == &other.elt_);
    return iota_it_ == other.iota_it_;
  }
  bool operator!=(const RepeatingIterator &other) const {
    return !operator==(other);
  }

  difference_type operator-(const RepeatingIterator &other) const {
    return iota_it_-other.iota_it_;
  }
  // auto dims_diff(const RepeatingIterator &other) const {
  //   return ::dims_diff(iota_it_, other.iota_it_);
  // }
};
template <typename Elt>
class RepeatingIterable {
public:
  using iterator = RepeatingIterator<Elt>;
  using const_iterator = RepeatingIterator<Elt>;

  const Elt elt_;
  const Iota<> iota_;

  explicit RepeatingIterable(const Elt &elt, const std::size_t &times) : elt_(elt), iota_(times) {
  }

  iterator begin() { return iterator(elt_, iota_.cbegin()); }
  iterator end() { return iterator(elt_, iota_.cend()); }
  const_iterator cbegin() const { return const_iterator(elt_, iota_.cbegin()); }
  const_iterator cend() const { return const_iterator(elt_, iota_.cend()); }
  const_iterator begin() const { return cbegin(); }
  const_iterator end() const { return cend(); }

  // xxx const iterators etc. for everything

  std::size_t size() const { return iota_.size(); }
};

// InferencesForTimesGenerator(curve, time_itbl) uses `curve` to infer the mean
// and standard deviation of the data at the times in `time_itbl`. If
// `time_itbl` coincides with the times in `curve`, then its output should be
// equivalent to `curve`.  This is like `approx` in R.
template<typename TimeIterable>
auto
InferencesForTimesGenerator(const ObservableCurve &curve, const TimeIterable &time_itbl) {
  auto ceiling_itbl = MakeCeilingIterable(time_itbl, curve,
            [](const std::pair<Time, Observable> &time_and_observable, Time time){
              return time_and_observable.first < time;
            });
  auto inference_itbl = MakeFMapIterable([&curve](const std::pair<proxy_for<iterator_for<const TimeIterable> >, iterator_for<const ObservableCurve> > &time_and_obsbl_ceil) {
      Time time = time_and_obsbl_ceil.first;
      const auto &obsbl_ceil = time_and_obsbl_ceil.second;
      const Observable &observable = InferWithCeiling(time, obsbl_ceil, curve);
      return observable;
    }, ceiling_itbl);
  return inference_itbl;
}

template<typename TimeIterable>
std::vector<Observable>
InferencesForTimesVector(const ObservableCurve &curve, const TimeIterable &time_itbl) {
  auto inference_itbl = InferencesForTimesGenerator(curve, time_itbl);
  std::vector<Observable> result(inference_itbl.begin(), inference_itbl.end());
  return result;
}

template<typename TimeIterable>
auto
IndexAtTimesGenerator(const Trajectory &trajectory, const TimeIterable &time_itbl) {
  auto ceiling_itbl = MakeCeilingIterable(time_itbl, trajectory,
            [](auto &&time_and_observation, Time time){
              return time_and_observation.first < time;
            });
  auto observation_itbl = MakeFMapIterable([](auto &&time_and_obs_ceil) {
      Time time = time_and_obs_ceil.first;
      const auto &obs_ceil = time_and_obs_ceil.second;
      Time ceil_time = obs_ceil->first;
      assert (time == ceil_time);
      if (time != ceil_time) throw std::logic_error("...");
      Observation obs = obs_ceil->second;
      return obs;
    }, ceiling_itbl);
  return observation_itbl;
}

template<std::size_t index>
struct GetFunctor {
  template<typename T>
  auto operator()(T &&t) const {
    using std::get;
    return get<index>(t);
  }
};

template<typename ShrinkageMap, typename ObservableIterable>
LogLikelihood
TrajectoryLogLikelihoodForAlignedObservables(const ShrinkageMap &shrinkage_map, const ObservableIterable &observable_itbl, const Trajectory &trajectory) {
  // xxx vs properly specifying type to select overload of std::get/ADL get
  auto time_itbl = MakeFMapIterable(GetFunctor<0>(), shrinkage_map);
  auto weight_itbl = MakeFMapIterable(GetFunctor<1>(), shrinkage_map);
  auto obs_itbl = IndexAtTimesGenerator(trajectory, time_itbl);
  auto log_lkhd_itbl = MakeFMapZipProductIterable(&Observable::log_likelihood, observable_itbl, obs_itbl);
  auto result = std::inner_product(weight_itbl.begin(), weight_itbl.end(), log_lkhd_itbl.begin(), LogLikelihood(0));
  return result;
}

template<typename ShrinkageMap>
LogLikelihood
TrajectoryLogLikelihood(const ShrinkageMap &shrinkage_map, const ObservableCurve &observable_curve, const Trajectory &trajectory) {
  // xxx vs properly specifying type to select overload of std::get/ADL get
  auto time_itbl = MakeFMapIterable(GetFunctor<0>(), shrinkage_map);
  auto observable_inference_itbl = InferencesForTimesGenerator(observable_curve, time_itbl);
  return TrajectoryLogLikelihoodForAlignedObservables(shrinkage_map, observable_inference_itbl, trajectory);
}

template<typename ShrinkageMap, typename TrajectoryIterable>
LogLikelihood
TrajectoryLogWeight(const ShrinkageMap &observed_past_shrinkage_map,
                    const ShrinkageMap &reasonable_future_shrinkage_map,
                    std::size_t n_future_neighbors,
                    const Trajectory &current_trajectory,
                    const TrajectoryIterable &historical_trajectories,
                    const ObservableCurve &observable_curve,
                    const Time bias_peaktime_mean,
                    const Time bias_peaktime_sd,
                    const Rp bias_peaktime_shrinkage,
                    const Mean bias_peakheight_mean,
                    const SD bias_peakheight_sd,
                    const Rp bias_peakheight_shrinkage
                    ) {
  if (( std::isnan(bias_peaktime_mean) ||  std::isnan(bias_peaktime_sd) ||  std::isnan(bias_peaktime_shrinkage)) &&
      (!std::isnan(bias_peaktime_mean) || !std::isnan(bias_peaktime_sd) || !std::isnan(bias_peaktime_shrinkage))) {
    throw std::invalid_argument("Either (a) all or (b) none of {bias_peaktime_mean, bias_peaktime_sd, bias_peaktime_shrinkage} must be NA.");
  }
  if (( std::isnan(bias_peakheight_mean) ||  std::isnan(bias_peakheight_sd) ||  std::isnan(bias_peakheight_shrinkage)) &&
      (!std::isnan(bias_peakheight_mean) || !std::isnan(bias_peakheight_sd) || !std::isnan(bias_peakheight_shrinkage))) {
    throw std::invalid_argument("Either (a) all or (b) none of {bias_peakheight_mean, bias_peakheight_sd, bias_peakheight_shrinkage} must be NA.");
  }
  LogLikelihood past_log_lkhd = TrajectoryLogLikelihood(observed_past_shrinkage_map, observable_curve, current_trajectory);
  auto future_times_itbl = MakeFMapIterable(GetFunctor<0>(), reasonable_future_shrinkage_map);
  std::vector<Observable> future_inferences = InferencesForTimesVector(observable_curve, future_times_itbl);
  std::vector<LogLikelihood> all_fut_wts;
  all_fut_wts.reserve(historical_trajectories.size());
  for (const auto &historical_trajectory : historical_trajectories) {
    all_fut_wts.push_back(TrajectoryLogLikelihoodForAlignedObservables(reasonable_future_shrinkage_map, future_inferences, historical_trajectory));
  }
  assert(n_future_neighbors <= all_fut_wts.size());
  auto nearest_fut_wt_end = all_fut_wts.begin()+static_cast<std::ptrdiff_t>(n_future_neighbors);
  std::partial_sort(all_fut_wts.begin(), nearest_fut_wt_end, all_fut_wts.end(),
                    std::greater<LogLikelihood>());
  LogLikelihood fut_log_weight = std::accumulate(all_fut_wts.begin(), nearest_fut_wt_end, LogLikelihood(0)) / Rp(n_future_neighbors);
  // Iterator pointing to the smoothed peak week and peak height:
  // fixme check whether smoothed peak is acceptable
  iterator_for<const ObservableCurve> peak_it = PeakIt(observable_curve);
  // Calculate posterior bias weight:
  LogLikelihood bias_weight = 0.0;
  if (!std::isnan(bias_peaktime_shrinkage)) {
    LogLikelihood bias_peaktime_weight = bias_peaktime_shrinkage * R::dnorm(peak_it->first, bias_peaktime_mean, bias_peaktime_sd, true);
    assert(!std::isnan(bias_peaktime_weight));
    bias_weight += bias_peaktime_weight;
  }
  if (!std::isnan(bias_peakheight_shrinkage)) {
    LogLikelihood bias_peakheight_weight = bias_peakheight_shrinkage * R::dnorm(peak_it->second.mean_, bias_peakheight_mean, bias_peakheight_sd, true);
    assert(!std::isnan(bias_peakheight_weight));
    bias_weight += bias_peakheight_weight;
  }
  return past_log_lkhd + fut_log_weight + bias_weight;
}

// xxx peak precomputation before yscale for peak
// xxx determine global+local yscale then perform

// xxx todo input validation checks for everything when inputted

// xxx rename most Cartesian stuff to "Outer"?

template<typename ShrinkageMap, typename ObservableCurveIterator>
auto LogWeightGenerator(const ShrinkageMap &observed_past_shrinkage_map,
                        const ShrinkageMap &reasonable_future_shrinkage_map,
                        std::size_t n_future_neighbors,
                        const Trajectory &current_trajectory,
                        const std::vector<Trajectory> &historical_trajectories,
                        const ObservableCurveIterator &obsbl_curve_it,
                        const Time bias_peaktime_mean,
                        const Time bias_peaktime_sd,
                        const Rp bias_peaktime_shrinkage,
                        const Mean bias_peakheight_mean,
                        const SD bias_peakheight_sd,
                        const Rp bias_peakheight_shrinkage
                        ) {
  return MakeFMapIterable([&, n_future_neighbors](const ObservableCurve &obsbl_curve){
    return TrajectoryLogWeight(observed_past_shrinkage_map,
                               reasonable_future_shrinkage_map,
                               n_future_neighbors,
                               current_trajectory,
                               historical_trajectories,
                               obsbl_curve,
                               bias_peaktime_mean,
                               bias_peaktime_sd,
                               bias_peaktime_shrinkage,
                               bias_peakheight_mean,
                               bias_peakheight_sd,
                               bias_peakheight_shrinkage
                               );
  }, obsbl_curve_it);
}

template<typename Products>
auto ObservableCurveGenerator(
const std::vector<ObservableCurve> &fit_observable_curves,
Mean y_scale_baseline,
const std::vector<std::ptrdiff_t> &curve_index_choices,
const std::vector<Time> &peak_time_choices,
const std::vector<Time> &x_shift_choices,
const std::vector<Time> &x_scale_choices,
const std::vector<SD> &sd_choices,
const std::vector<Rp> &sd_scale_choices,
const std::vector<Mean> &peak_height_choices,
const std::vector<Rp> &y_scale_choices
) {
  // todo std::move
  auto &&fit_curve_ptr_choices =
      MakeShallowCacheIterable(MakeFMapIterable(
          [&fit_observable_curves](std::ptrdiff_t curve_index) {
            return &fit_observable_curves[static_cast<std::size_t>(curve_index)];
          },
          curve_index_choices));
  auto &&fit_curve_xformers =
      MakeShallowCacheIterable(MakeFMapIterable(
          [](const ObservableCurve *fit_curve_ptr_choice) {
            return MakeXTransformer(*fit_curve_ptr_choice);
          },
          fit_curve_ptr_choices));
  auto &&global_x_shifted_xformers = MakeShallowCacheIterable(
      Products::MakeFMap(&XTransformer::XShiftedForPeakTime, fit_curve_xformers,
                         peak_time_choices));
  auto &&local_x_shifted_xformers = MakeShallowCacheIterable(
      Products::MakeFMap(&XTransformer::XShiftedRight,
                         global_x_shifted_xformers, x_shift_choices));
  auto &&x_scaled_xformers = MakeShallowCacheIterable(
      Products::MakeFMap(&XTransformer::XScaledWidenAroundPeakTime,
                         local_x_shifted_xformers, x_scale_choices));
  auto &&x_scaled_curves = MakeShallowCacheIterable(MakeFMapIterable(&XTransformer::TransformedCurve, x_scaled_xformers));
  auto &&sd_assigned_curves = MakeShallowCacheIterable(Products::MakeFMap(
                  &ObservableCurveWithSD, x_scaled_curves, sd_choices));
  auto &&sd_scaled_curves = MakeShallowCacheIterable(Products::MakeFMap(
                  &ObservableCurveWithSDScale, sd_assigned_curves, sd_scale_choices));
  auto &&global_y_scaled_curves = MakeShallowCacheIterable(Products::MakeFMap(PeakHeightScalerAboveBaseline(y_scale_baseline), sd_scaled_curves, peak_height_choices));
  auto &&local_y_scaled_curves = MakeShallowCacheIterable(Products::MakeFMap(YScalerAboveBaseline(y_scale_baseline), global_y_scaled_curves, y_scale_choices));
  return local_y_scaled_curves;
}

template <typename Products>
auto ExampleTransformedCurves(const std::vector<ObservableCurve> &fit_observable_curves//,
            // const std::vector<std::ptrdiff_t> &curve_index_choices,
            // const std::vector<Time> &peak_time_choices
            ) {
  Mean y_scale_baseline = 1.2;

  std::vector<std::ptrdiff_t> curve_index_choices = {{0, 1}};
  std::vector<Time> peak_time_choices = {{10, NAN}};
  // std::vector<Time> peak_time_choices = {{7, 8, 9, 10, 11, 12, 13, 14, 15, 16, NAN}};
  std::vector<Time> x_shift_choices = {{-1, 0}};
  // std::vector<Time> x_shift_choices = {{-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5}};
  std::vector<Rp> x_scale_choices = {{1, 2}};
  // std::vector<Rp> x_scale_choices = {{1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2}};
  std::vector<SD> sd_choices = {{0.5, NAN}};
  std::vector<Rp> sd_scale_choices = {{1, 2}};
  std::vector<Mean> peak_height_choices = {{NAN, 10}};
  std::vector<Rp> y_scale_choices = {{0.9,1.1}};

  return ObservableCurveGenerator<Products>(
          fit_observable_curves,
          y_scale_baseline,
          curve_index_choices,
          peak_time_choices,
          x_shift_choices,
          x_scale_choices,
          sd_choices,
          sd_scale_choices,
          peak_height_choices,
          y_scale_choices
          );
}

// xxx enable_if on product templates / constructors

// xxx reordering of transformations and optimization of zip products for repeated discrete choices

// xxx should not shallow cache zip results; have ShallowCacheCartesianProducts
// instead

// int main() {
//   // maino();
//   std::vector<MeanCurve> fit_mean_curves = {
//     {{0, 4}, {1, 7}, {2, 5}, {3, 3}, {4, 2}},
//     {{0, 1}, {1, 2}, {2, 4}, {3, 3}, {4, 3}, {5, 1.8}}};
//   std::vector<Time> points = {{-1,-0.5,0,0.5,1,1.5,1.8,1.9,2,2,2,2.5,3,4,10}};
//   std::vector<Time> levels = {{0,1,2,2,2,2.6}};
//   const Trajectory trajectory = {{{-1,1}, {0,2}, {1,4}, {2,5}, {3,6}, {4,3}, {5,2}}};
//   // const std::vector<std::pair<Time, Rp>> shrinkage_map = {{{0,1},{1,1},{2,1},{3,1},{4,1},{5,1}}};
//   const std::vector<std::pair<Time, Rp>> shrinkage_map = {{{0,1},{1,1},{2,1},{3,1},{4,1}}};
//   // auto time_itbl = MakeFMapIterable(GetFunctor<0>(), shrinkage_map);
//   // auto weight_itbl = MakeFMapIterable(GetFunctor<1>(), shrinkage_map);
//   // auto obs_itbl = IndexAtTimesGenerator(trajectory, time_itbl);
//   // std::cout << "TRAJECTORY AT SHRINKAGE MAP TIMES:\n";
//   // for (auto &&asdf : obs_itbl) {
//   //   std::cout << asdf << "\n";
//   // }
//   // auto shrinkage_map = MakeZipProductIterable(Iota<Time>(5), RepeatingIterable<Rp>(1, 5));
//   // std::vector<std::pair<Time, Rp>> shrinkage_map = {{{0,1},{1,1},{2,1},{3,1},{4,1},{5,1}}};
//   // auto shrinkage_map = MakeZipProductIterable(Iota<Time>(5), RepeatingIterable<Rp>(1, 5));
//   // std::cout << "WEIGHT MAP:" << std::endl;
//   // for (auto &&shrinkage_mapping : shrinkage_map) {
//   //   std::cout << std::get<0>(shrinkage_mapping) << ": " << std::get<1>(shrinkage_mapping) << std::endl;
//   // }
//   std::vector<SD> fit_sds = {{0.1,0.8}};
//   auto fit_observable_curves_gen = MakeFMapZipProductIterable(&MakeObservableCurve, fit_mean_curves, fit_sds);
//   const std::vector<ObservableCurve> fit_observable_curves(fit_observable_curves_gen.begin(), fit_observable_curves_gen.end());
//   // ObservableCurve curve = fit_observable_curves[0];
//   // auto obsbl_itbl = InferencesForTimesGenerator(curve, time_itbl);
//   // auto log_lkhd_itbl = MakeFMapZipProductIterable(&Observable::log_likelihood, obsbl_itbl, obs_itbl);
//   // std::cout << "ll's" << std::endl;
//   // for (LogLikelihood ll : log_lkhd_itbl) {
//   //   std::cout << ll << "\n";
//   // }
//   // std::cout << "INFER" << std::endl;
//   // for (const auto &asdf : MakeTrajectoryLogLikelihoodsIterator(trajectory, fit_observable_curves[0])) {
//   //   std::cout << asdf << std::endl;
//   // }
//   // const auto &time_itbl_pre = MakeFMapIterable(GetFunctor<0>(), shrinkage_map);
//   // auto ceiling_itbl = MakeCeilingIterable(time_itbl_pre, fit_observable_curves[0],
//   //            [](const std::pair<Time, Observable> &time_and_observable, Time time){
//   //              return time_and_observable.first < time;
//   //            });
//   // for (auto &&asdf : ceiling_itbl) {
//   //   std::cout << asdf.first << std::endl;
//   // }
//   // auto curve = fit_observable_curves[0];
//   // using TimeIterable = typename std::remove_reference<decltype(time_itbl_pre)>::type;
//   // using ObservableCurve = typename std::remove_reference<decltype(curve)>::type;
//   // const auto &inference_itbl = MakeFMapIterable([&curve](const std::pair<
//   //              proxy_for<iterator_for<const TimeIterable> >,
//   //              iterator_for<const ObservableCurve> >
//   //              &time_and_obsbl_ceil) {
//   // // auto inference_itbl = MakeFMapIterable([&curve](const std::pair<
//   // //               proxy_for<iterator_for<const TimeIterable> >,
//   // //               iterator_for<const ObservableCurve> >
//   // //               &time_and_obsbl_ceil) {
//   // // const auto &inference_itbl = MakeFMapIterable([&curve](const auto &time_and_obsbl_ceil) {
//   //     Time time = time_and_obsbl_ceil.first;
//   //     const auto &obsbl_ceil = time_and_obsbl_ceil.second;
//   //     const Observable &observable = InferWithCeiling(time, obsbl_ceil, curve);
//   //     return observable;
//   //   }, ceiling_itbl);
//   // for (auto &&asdf : inference_itbl) {
//   //   std::cout << asdf << std::endl;
//   // }
//   // const auto &observable_inference_itbl = InferencesForTimesGenerator(fit_observable_curves[0], time_itbl_pre);
//   // auto &&observable_inference_itbl = InferencesForTimesGenerator(fit_observable_curves[0], time_itbl_pre);
//   // for (const Observable &asdf : observable_inference_itbl) {
//   //   std::cout << asdf << std::endl;
//   // }
//   // auto time_itbl = MakeFMapIterable(GetFunctor<0>(), shrinkage_map);
//   // auto weight_itbl = MakeFMapIterable(GetFunctor<1>(), shrinkage_map);
//   // auto obs_itbl = IndexAtTimesGenerator(trajectory, time_itbl);
//   // auto log_lkhd_itbl = MakeFMapZipProductIterable(&Observable::log_likelihood, observable_inference_itbl, obs_itbl);
//   // auto ll = std::inner_product(weight_itbl.begin(), weight_itbl.end(), log_lkhd_itbl.begin(), LogLikelihood(0));
//   // std::cout << "LL: " << ll << std::endl;





//   std::cout << "LL: " << TrajectoryLogLikelihood(shrinkage_map, fit_observable_curves[0], trajectory) << std::endl;
//   std::cout << "LW: " << TrajectoryLogWeight(shrinkage_map, shrinkage_map, 2,
//                  trajectory, fit_mean_curves,
//                  fit_observable_curves[0]
//                ) << std::endl;

//   // std::vector<std::ptrdiff_t> curve_index_choices = {{0, 1}};
//   // std::vector<Time> peak_time_choices = {{10, NAN}};
//   std::cout << "CART" << std::endl;
//   for (const auto &xformed_curve :
//          ExampleTransformedCurves<CartesianProducts>(fit_observable_curves/*, curve_index_choices*//*, peak_time_choices*/)) {
//     std::cout << XTransformer(xformed_curve) << std::endl;
//   }
//   double x = 0;
//   const auto &xformed_curves = ExampleTransformedCurves<CartesianProducts>(fit_observable_curves/*, curve_index_choices*//*, peak_time_choices*/);
//   std::cout << "CARTinitd" << std::endl;
//   for (const auto &xformed_curve : xformed_curves) {
//     // Ignore(xformed_curve);
//     // std::cout << XTransformer(xformed_curve) << std::endl;
//     // ++x;
//     // x += xformed_curve.size();
//     x += xformed_curve.begin()->second.mean_;
//   }
//   std::cout << xformed_curves.size() << " " << x << std::endl;

//   for (LogLikelihood log_weight : LogWeightGenerator(
//           shrinkage_map,
//           shrinkage_map,
//           2,
//           trajectory,
//           fit_mean_curves,
//           xformed_curves
//                    )) {
//     std::cout << log_weight << "\n";
//   }

//   std::cout << "ZIP" << std::endl;
//   for (const auto &xformed_curve :
//          ExampleTransformedCurves<ZipProducts>(fit_observable_curves/*, curve_index_choices*//*, peak_time_choices*/)) {
//     std::cout << XTransformer(xformed_curve) << std::endl;
//   }
//   // Iota<Time> times(5);
//   // for (auto time : times) {
//   //   std::cout << time << std::endl;
//   // }
//   // std::cout << times.end() - times.begin() << std::endl;
//   // for (auto time : times) {
//   //   std::cout << time << std::endl;
//   // }
//   // std::cout << times.end() - times.begin() << std::endl;
//   // auto time_it = times.begin();
//   // ++time_it;
//   // std::cout << time_it - times.begin() << std::endl;
//   // std::cout << times.begin() - time_it << std::endl;
//   // std::cout << "EXIT MAIN" << std::endl;

//   OUTELT(Iota<>(5).ndims());
//   OUTELT(Iota<>::ndims());
//   OUTELT(ndims(Iota<>(5)));
//   // OUTELT(ndims<Iota<>>());
//   OUTELT(t_ndims<Iota<>>);
//   OUTELT(NDims<Iota<>>::value);
//   OUTELT(std::get<0>(Iota<>(5).dims()));
//   OUTELT(std::get<0>(dims(Iota<>(5))));
//   OUTITBL(dims(Iota<>(5)));
//   // OUTITBL(dims(xformed_curves));
//   OUTITBL(dims(fit_observable_curves_gen));
//   std::vector<int> a{{1,2}};
//   std::vector<double> b{{0.1,0.2,0.3}};
//   auto &&ab = MakeCartesianProductIterable(a,b);
//   OUTELT(ndims(ab));
//   OUTITBL(dims(ab));
//   OUTITBL(dims(xformed_curves));

//   return 0;
// }

// todo iterable references / uniqueptr's / etc. --- forwarding... prevent copying into other iterables
// todo "scalar" wrappers
// xxx FMapZip --> just FMap

// Doesn't seem to have an effect:
// [[Rcpp::plugins(cpp14)]]

RcppExport SEXP ExampleZipCurves() {
  BEGIN_RCPP
    // maino();
    std::vector<MeanCurve> fit_mean_curves = {
        {{0, 4}, {1, 7}, {2, 5}, {3, 3}, {4, 2}},
        {{0, 1}, {1, 2}, {2, 4}, {3, 3}, {4, 3}, {5, 1.8}}};
    std::vector<SD> fit_sds = {{0.1,0.8}};
    auto fit_observable_curves_gen = MakeFMapZipProductIterable(&MakeObservableCurve, fit_mean_curves, fit_sds);
    std::vector<ObservableCurve> fit_observable_curves(fit_observable_curves_gen.begin(), fit_observable_curves_gen.end());
    auto xformed_observable_curves = ExampleTransformedCurves<ZipProducts>(fit_observable_curves);
    auto xformed_mean_maps_gen = MakeFMapIterable([](const ObservableCurve &observable_curve){
  auto mean_curve_mappings = MakeFMapIterable([](const proxy_for<iterator_for<decltype(observable_curve)> > &mapping) {
      return std::make_pair(mapping.first,mapping.second.mean_);
    }, observable_curve);
  return std::map<Time, Mean>(mean_curve_mappings.begin(), mean_curve_mappings.end());
      }, xformed_observable_curves);
    std::vector<std::map<Time, Mean> > xformed_mean_maps(xformed_mean_maps_gen.begin(), xformed_mean_maps_gen.end());
    return Rcpp::wrap(xformed_mean_maps);
  END_RCPP
}

RcppExport SEXP ExampleCartesianCurves() {
  BEGIN_RCPP
    // maino();
    std::vector<MeanCurve> fit_mean_curves = {
        {{0, 4}, {1, 7}, {2, 5}, {3, 3}, {4, 2}},
        {{0, 1}, {1, 2}, {2, 4}, {3, 3}, {4, 3}, {5, 1.8}}};
    std::vector<SD> fit_sds = {{0.1,0.8}};
    auto fit_observable_curves_gen = MakeFMapZipProductIterable(&MakeObservableCurve, fit_mean_curves, fit_sds);
    std::vector<ObservableCurve> fit_observable_curves(fit_observable_curves_gen.begin(), fit_observable_curves_gen.end());
    auto xformed_observable_curves = ExampleTransformedCurves<CartesianProducts>(fit_observable_curves);
    auto xformed_mean_maps_gen = MakeFMapIterable([](const ObservableCurve &observable_curve){
  auto mean_curve_mappings = MakeFMapIterable([](const proxy_for<iterator_for<decltype(observable_curve)> > &mapping) {
      return std::make_pair(mapping.first,mapping.second.mean_);
    }, observable_curve);
  return std::map<Time, Mean>(mean_curve_mappings.begin(), mean_curve_mappings.end());
      }, xformed_observable_curves);
    std::vector<std::map<Time, Mean> > xformed_mean_maps(xformed_mean_maps_gen.begin(), xformed_mean_maps_gen.end());
    return Rcpp::wrap(xformed_mean_maps);
  END_RCPP
}

Trajectory TrajectoryFromDataFrame(const Rcpp::DataFrame &trajectory_df) {
  // todo auto &&, move, ...
  auto trajectory_times = Rcpp::as<std::vector<Time>>(trajectory_df["time"]);
  auto trajectory_observations = Rcpp::as<std::vector<Observation>>(trajectory_df["observation"]);
  auto trajectory_gen = MakeFMapZipProductIterable(&std::make_pair<const Time&, const Observation&>, trajectory_times, trajectory_observations);
  auto trajectory = Trajectory(trajectory_gen.begin(), trajectory_gen.end());
  return trajectory;
}

std::vector<std::pair<Time, Rp>> ShrinkageMapFromDataFrame(const Rcpp::DataFrame &shrinkage_map_df) {
  auto times = Rcpp::as<std::vector<Time>>(shrinkage_map_df["time"]);
  auto weights = Rcpp::as<std::vector<Rp>>(shrinkage_map_df["weight"]);
  auto gen = MakeFMapZipProductIterable(&std::make_pair<const Time&, const Rp&>, times, weights);
  std::vector<std::pair<Time, Rp>> result(gen.begin(), gen.end());
  return result;
}

// todo pack arguments into a list, convert to POD struct and pass around that way.

// [[Rcpp::export]]
std::vector<Rcpp::DataFrame> CartesianProductCurves
(
 Rcpp::List fit_obj_rcpp,
 Mean y_scale_baseline,
 std::vector<std::ptrdiff_t> curve_index_choices,
 std::vector<Time> peak_time_choices,
 std::vector<Time> x_shift_choices,
 std::vector<Time> x_scale_choices,
 std::vector<SD> sd_choices,
 std::vector<Rp> sd_scale_choices,
 std::vector<Mean> peak_height_choices,
 std::vector<Rp> y_scale_choices) {
    // ProfilerStart("./myprofile.log");
      // xxx prevent copying using proxies to Rcpp versions
  std::vector<ObservableCurve> fit_observable_curves;
  for (const SEXP &curve_sexp : fit_obj_rcpp) {
    Rcpp::List curve_rcpp(curve_sexp);
    const auto type = Rcpp::as<std::string>(curve_rcpp["type"]);
    if (type != "Gaussian")
      throw std::invalid_argument(
          "Each fit curve in fit.obj must have $type==\"Gaussian\"");
    const auto f = Rcpp::as<std::vector<Mean>>(curve_rcpp["f"]);
    const auto tau = Rcpp::as<SD>(curve_rcpp["tau"]);
    const Iota<Time> curve_times(f.size());
    auto mean_curve_gen = MakeFMapZipProductIterable(&std::make_pair<Time, const Mean &>, curve_times, f);
    MeanCurve mean_curve(mean_curve_gen.begin(), mean_curve_gen.end());
    ObservableCurve observable_curve = MakeObservableCurve(mean_curve, tau);
    fit_observable_curves.push_back(observable_curve);
  }

  // Rcpp::Rcout << "make result" << std::endl;
  const auto &result_gen = ObservableCurveGenerator<CartesianProducts>(
               // new_dat,
               // dat,
               fit_observable_curves,
               y_scale_baseline,
               // observed_past_shrinkage_map,
               // reasonable_future_shrinkage_map,
               // n_future_neighbors,
               curve_index_choices,
               peak_time_choices,
               x_shift_choices,
               x_scale_choices,
               sd_choices,
               sd_scale_choices,
               peak_height_choices,
               y_scale_choices
               );
  // std::vector<ObservableCurve> result_vec(result.begin(), result.end());
  // ProfilerStop();
  // return Rcpp::wrap(result_vec);
  // Rcpp::List result_rcpp(result.size());
  // Rcpp::IntegerVector result_dims_vec(ndims(result));
  // const auto &result_dims_ray = dims(result);
  // std::copy(result_dims_ray.begin(), result_dims_ray.end(), result_dims_vec);
  // fixme must reverse, and doesn't work anyway
  // auto ray_it = result_dims_ray.begin(), ray_end = result_dims_ray.end();
  // for (auto vec_it = result_dims_vec.begin(); ray_it != ray_end; ++ray_it, ++vec_it) {
  //   *vec_it = *ray_it;
  // }
  // Rcpp::Dimension result_dims(result_dims_vec);
  // Rcpp::List result_rcpp(result_dims);
  std::vector<Rcpp::DataFrame> result_vec;
  result_vec.reserve(result_gen.size());
  // Rcpp::Rcout << "ncurves: " << result_gen.size() << std::endl;
  // for (auto xformed_curve : result) {
  // std::size_t ind=0;
  // for (auto it = result.begin(), end = result.end(); it != end; ++it, ++ind) {
  // for (auto it = result_gen.cbegin(), end = result_gen.cend(); it != end; ++it, ++ind) {
  for (const auto &xformed_curve : result_gen) {
    // fixme why does this cause conversion errors?  both should be const iterators...
    // Rcpp::Rcout << ind << ": " << (end-it) << "= -" << (it-end) << " remaining...\n";
    std::vector<Time> times;
    std::vector<Mean> means;
    std::vector<SD> sds;
    times.reserve(xformed_curve.size());
    means.reserve(xformed_curve.size());
    sds.reserve(xformed_curve.size());
    for (auto time_and_observable : xformed_curve) {
      times.push_back(time_and_observable.first);
      means.push_back(time_and_observable.second.mean_);
      sds.push_back(time_and_observable.second.sd_);
      // Rcpp::Rcout << time_and_observable.first << ": " << time_and_observable.second << "\n";
    }
    // Rcpp::Rcout << std::endl;
    Rcpp::DataFrame curve_rcpp = Rcpp::DataFrame::create(
               Rcpp::Named("time")=Rcpp::wrap(times),
               Rcpp::Named("mean")=Rcpp::wrap(means),
               Rcpp::Named("sd")=Rcpp::wrap(sds)
               );
    result_vec.push_back(curve_rcpp);
  }
  return result_vec;
}

// [[Rcpp::export]]
Rcpp::NumericVector CartesianProductLogWeights
(
 Rcpp::DataFrame new_dat_df,
 Rcpp::List dat_rcpp,
 Rcpp::DataFrame observed_past_shrinkage_map_df,
 Rcpp::DataFrame reasonable_future_shrinkage_map_df,
 std::size_t n_future_neighbors,
 Rcpp::List fit_obj_rcpp,
 Mean y_scale_baseline,
 std::vector<std::ptrdiff_t> curve_index_choices,
 std::vector<Time> peak_time_choices,
 std::vector<Time> x_shift_choices,
 std::vector<Time> x_scale_choices,
 std::vector<SD> sd_choices,
 std::vector<Rp> sd_scale_choices,
 std::vector<Mean> peak_height_choices,
 std::vector<Rp> y_scale_choices,
 Time bias_peaktime_mean,
 Time bias_peaktime_sd,
 Rp bias_peaktime_shrinkage,
 Mean bias_peakheight_mean,
 SD bias_peakheight_sd,
 Rp bias_peakheight_shrinkage
 ) {
    // ProfilerStart("./myprofile.log");
      // xxx prevent copying using proxies to Rcpp versions
  Trajectory current_trajectory = TrajectoryFromDataFrame(new_dat_df);
  std::vector<Trajectory> historical_trajectories;
  historical_trajectories.reserve(dat_rcpp.length());
  for (const SEXP &historical_trajectory_sexp : dat_rcpp) {
    historical_trajectories.push_back(TrajectoryFromDataFrame(Rcpp::DataFrame(historical_trajectory_sexp)));
  }
  // Rcpp::Rcout << "observed_past_shrinkage_map" << std::endl;
  std::vector<std::pair<Time,Rp>> observed_past_shrinkage_map = ShrinkageMapFromDataFrame(observed_past_shrinkage_map_df);
  // Rcpp::Rcout << "reasonable_future_shrinkage_map" << std::endl;
  std::vector<std::pair<Time,Rp>> reasonable_future_shrinkage_map = ShrinkageMapFromDataFrame(reasonable_future_shrinkage_map_df);
  // Rcpp::Rcout << "fit.obj" << std::endl;
  std::vector<ObservableCurve> fit_observable_curves;
  for (const SEXP &curve_sexp : fit_obj_rcpp) {
    Rcpp::List curve_rcpp(curve_sexp);
    const auto type = Rcpp::as<std::string>(curve_rcpp["type"]);
    if (type != "Gaussian")
      throw std::invalid_argument(
          "Each fit curve in fit.obj must have $type==\"Gaussian\"");
    const auto f = Rcpp::as<std::vector<Mean>>(curve_rcpp["f"]);
    const auto tau = Rcpp::as<SD>(curve_rcpp["tau"]);
    const Iota<Time> curve_times(f.size());
    auto mean_curve_gen = MakeFMapZipProductIterable(&std::make_pair<Time, const Mean &>, curve_times, f);
    MeanCurve mean_curve(mean_curve_gen.begin(), mean_curve_gen.end());
    ObservableCurve observable_curve = MakeObservableCurve(mean_curve, tau);
    fit_observable_curves.push_back(observable_curve);
  }

  // Rcpp::Rcout << "CART make result" << std::endl;
  auto &&xformed_curves_gen = ObservableCurveGenerator<CartesianProducts>(
               // new_dat,
               // dat,
               fit_observable_curves,
               y_scale_baseline,
               // observed_past_shrinkage_map,
               // reasonable_future_shrinkage_map,
               // n_future_neighbors,
               curve_index_choices,
               peak_time_choices,
               x_shift_choices,
               x_scale_choices,
               sd_choices,
               sd_scale_choices,
               peak_height_choices,
               y_scale_choices
               );
  auto &&log_weights_gen = LogWeightGenerator(observed_past_shrinkage_map,
                                              reasonable_future_shrinkage_map,
                                              n_future_neighbors,
                                              current_trajectory,
                                              historical_trajectories,
                                              std::move(xformed_curves_gen),
                                              bias_peaktime_mean,
                                              bias_peaktime_sd,
                                              bias_peaktime_shrinkage,
                                              bias_peakheight_mean,
                                              bias_peakheight_sd,
                                              bias_peakheight_shrinkage
                                              );
  // Rcpp::Rcout << "CART make result gen initd" << std::endl;
  // std::vector<ObservableCurve> result_vec(result.begin(), result.end());
  // ProfilerStop();
  // return Rcpp::wrap(result_vec);
  // Rcpp::List result_rcpp(result.size());
  // Rcpp::IntegerVector result_dims_vec(ndims(result));
  // const auto &result_dims_ray = dims(result);
  // std::copy(result_dims_ray.begin(), result_dims_ray.end(), result_dims_vec);
  // fixme must reverse, and doesn't work anyway
  // auto ray_it = result_dims_ray.begin(), ray_end = result_dims_ray.end();
  // for (auto vec_it = result_dims_vec.begin(); ray_it != ray_end; ++ray_it, ++vec_it) {
  //   *vec_it = *ray_it;
  // }
  // Rcpp::Dimension result_dims(result_dims_vec);
  // Rcpp::List result_rcpp(result_dims);
  // xxx vs init with size and set.  check distance operators...
  Rcpp::NumericVector result_rcpp(log_weights_gen.begin(), log_weights_gen.end());
  return result_rcpp;
}

// RcppExport SEXP ZipProductCurvesAndLogWeightsRcpp
// (
//  SEXP new_dat_sexp,
//  SEXP dat_sexp,
//  SEXP observed_past_shrinkage_map_sexp,
//  SEXP reasonable_future_shrinkage_map_sexp,
//  SEXP n_future_neighbors_sexp,
//  SEXP fit_obj_sexp,
//  SEXP y_scale_baseline_sexp,
//  SEXP curve_index_choices_sexp,
//  SEXP peak_time_choices_sexp,
//  SEXP x_shift_choices_sexp,
//  SEXP x_scale_choices_sexp,
//  SEXP sd_choices_sexp,
//  SEXP sd_scale_choices_sexp,
//  SEXP peak_height_choices_sexp,
//  SEXP y_scale_choices_sexp) {
//   BEGIN_RCPP  // forward exceptions to R
//     // ProfilerStart("./myprofile.log");
//       // xxx prevent copying using proxies to Rcpp versions
//   Rcpp::Rcout << "new.dat" << std::endl;
//   Trajectory current_trajectory = TrajectoryFromSEXP(new_dat_sexp);
//   Rcpp::Rcout << "dat" << std::endl;
//   auto dat_rcpp = Rcpp::List(dat_sexp);
//   std::vector<Trajectory> historical_trajectories;
//   historical_trajectories.reserve(dat_rcpp.length());
//   for (const SEXP &historical_trajectory_sexp : dat_rcpp) {
//     historical_trajectories.push_back(TrajectoryFromSEXP(historical_trajectory_sexp));
//   }
//   Rcpp::Rcout << "observed_past_shrinkage_map" << std::endl;
//   std::vector<std::pair<Time,Rp>> observed_past_shrinkage_map = ShrinkageMapFromSEXP(observed_past_shrinkage_map_sexp);
//   Rcpp::Rcout << "reasonable_future_shrinkage_map" << std::endl;
//   std::vector<std::pair<Time,Rp>> reasonable_future_shrinkage_map = ShrinkageMapFromSEXP(reasonable_future_shrinkage_map_sexp);
//   auto n_future_neighbors = Rcpp::as<std::size_t>(n_future_neighbors_sexp);
//   Rcpp::Rcout << "fit.obj" << std::endl;
//   Rcpp::List fit_obj_rcpp(fit_obj_sexp);
//   std::vector<ObservableCurve> fit_observable_curves;
//   for (const SEXP &curve_sexp : fit_obj_rcpp) {
//     Rcpp::List curve_rcpp(curve_sexp);
//     const auto type = Rcpp::as<std::string>(curve_rcpp["type"]);
//     if (type != "Gaussian")
//       throw std::invalid_argument(
//           "Each fit curve in fit.obj must have $type==\"Gaussian\"");
//     const auto f = Rcpp::as<std::vector<Mean>>(curve_rcpp["f"]);
//     const auto tau = Rcpp::as<SD>(curve_rcpp["tau"]);
//     const Iota<Time> curve_times(f.size());
//     auto mean_curve_gen = MakeFMapZipProductIterable(&std::make_pair<Time, const Mean &>, curve_times, f);
//     MeanCurve mean_curve(mean_curve_gen.begin(), mean_curve_gen.end());
//     ObservableCurve observable_curve = MakeObservableCurve(mean_curve, tau);
//     fit_observable_curves.push_back(observable_curve);
//   }

//   auto y_scale_baseline = Rcpp::as<Mean>(y_scale_baseline_sexp);

//   Rcpp::Rcout << "choices" << std::endl;
//   auto curve_index_choices = Rcpp::as<std::vector<std::ptrdiff_t>>(curve_index_choices_sexp);
//   auto peak_time_choices = Rcpp::as<std::vector<Time>>(peak_time_choices_sexp);
//   auto x_shift_choices = Rcpp::as<std::vector<Time>>(x_shift_choices_sexp);
//   auto x_scale_choices = Rcpp::as<std::vector<Time>>(x_scale_choices_sexp);
//   auto sd_choices = Rcpp::as<std::vector<SD>>(sd_choices_sexp);
//   auto sd_scale_choices = Rcpp::as<std::vector<Rp>>(sd_scale_choices_sexp);
//   auto peak_height_choices = Rcpp::as<std::vector<Mean>>(peak_height_choices_sexp);
//   auto y_scale_choices = Rcpp::as<std::vector<Rp>>(y_scale_choices_sexp);

//   Rcpp::Rcout << "ZIP make result" << std::endl;
//   auto &&xformed_curves_gen = ObservableCurveGenerator<ZipProducts>(
// 							 // new_dat,
// 							 // dat,
// 							 fit_observable_curves,
// 							 y_scale_baseline,
// 							 // observed_past_shrinkage_map,
// 							 // reasonable_future_shrinkage_map,
// 							 // n_future_neighbors,
// 							 curve_index_choices,
// 							 peak_time_choices,
// 							 x_shift_choices,
// 							 x_scale_choices,
// 							 sd_choices,
// 							 sd_scale_choices,
// 							 peak_height_choices,
// 							 y_scale_choices
// 							 );
//   Rcpp::Rcout << "make result gen initd" << std::endl;
//   // const auto &log_weights_gen = LogWeightGenerator(
//   // 						   observed_past_shrinkage_map,
//   // 						   reasonable_future_shrinkage_map,
//   // 						   n_future_neighbors,
//   // 						   current_trajectory,
//   // 						   historical_trajectories,
//   // 						   xformed_curves_gen);
//   // Rcpp::List curves_rcpp(xformed_curves_gen.size());
//   // Rcpp::NumericVector log_weights_rcpp(xformed_curves_gen.size());
//   std::vector<Rcpp::DataFrame> curves_vec;
//   std::vector<LogLikelihood> log_weights_vec;
//   curves_vec.reserve(xformed_curves_gen.size());
//   log_weights_vec.reserve(xformed_curves_gen.size());
//   std::size_t count = 0;
//   for (const auto &xformed_curve : xformed_curves_gen) {
//     if (count % 200 == 0) {
//       Rcpp::Rcout << count << std::endl;
//     }
//     ++count;
//     std::vector<Time> times;
//     std::vector<Mean> means;
//     std::vector<SD> sds;
//     times.reserve(xformed_curve.size());
//     means.reserve(xformed_curve.size());
//     sds.reserve(xformed_curve.size());
//     for (auto time_and_observable : xformed_curve) {
//       times.push_back(time_and_observable.first);
//       means.push_back(time_and_observable.second.mean_);
//       sds.push_back(time_and_observable.second.sd_);
//       // Rcpp::Rcout << time_and_observable.first << ": " << time_and_observable.second << "\n";
//     }
//     // Rcpp::Rcout << std::endl;
//     Rcpp::DataFrame curve_rcpp = Rcpp::DataFrame::create(
// 							 Rcpp::Named("time")=Rcpp::wrap(std::move(times)),
// 							 Rcpp::Named("mean")=Rcpp::wrap(std::move(means)),
// 							 Rcpp::Named("sd")=Rcpp::wrap(std::move(sds))
// 							 );
//     curves_vec.push_back(curve_rcpp);
//     log_weights_vec.push_back(TrajectoryLogWeight(
// 						  observed_past_shrinkage_map,
// 						  reasonable_future_shrinkage_map,
// 						  n_future_neighbors,
// 						  current_trajectory,
// 						  historical_trajectories,
// 						  xformed_curve));
//   }
//   Rcpp::List result_rcpp({{Rcpp::wrap(curves_vec), Rcpp::wrap(log_weights_vec)}});
//   return result_rcpp;
//   END_RCPP  // forward exceptions to R
// }

// [[Rcpp::export]]
SEXP ZipProductCurvesAndLogWeightsp
(
 std::vector<Time> output_times,
 Rcpp::DataFrame new_dat_df,
 Rcpp::List dat_rcpp,
 Rcpp::DataFrame observed_past_shrinkage_map_df,
 Rcpp::DataFrame reasonable_future_shrinkage_map_df,
 std::size_t n_future_neighbors,
 Rcpp::List fit_obj_rcpp,
 Mean y_scale_baseline,
 std::vector<std::ptrdiff_t> curve_index_choices,
 std::vector<Time> peak_time_choices,
 std::vector<Time> x_shift_choices,
 std::vector<Time> x_scale_choices,
 std::vector<SD> sd_choices,
 std::vector<Rp> sd_scale_choices,
 std::vector<Mean> peak_height_choices,
 std::vector<Rp> y_scale_choices,
 Time bias_peaktime_mean,
 Time bias_peaktime_sd,
 Rp bias_peaktime_shrinkage,
 Mean bias_peakheight_mean,
 SD bias_peakheight_sd,
 Rp bias_peakheight_shrinkage
 ) {
    // ProfilerStart("./myprofile.log");
      // xxx prevent copying using proxies to Rcpp versions
  Trajectory current_trajectory = TrajectoryFromDataFrame(new_dat_df);
  std::vector<Trajectory> historical_trajectories;
  historical_trajectories.reserve(dat_rcpp.length());
  for (const SEXP &historical_trajectory_sexp : dat_rcpp) {
    historical_trajectories.push_back(TrajectoryFromDataFrame(Rcpp::DataFrame(historical_trajectory_sexp)));
  }
  std::vector<std::pair<Time,Rp>> observed_past_shrinkage_map = ShrinkageMapFromDataFrame(observed_past_shrinkage_map_df);
  std::vector<std::pair<Time,Rp>> reasonable_future_shrinkage_map = ShrinkageMapFromDataFrame(reasonable_future_shrinkage_map_df);
  std::vector<ObservableCurve> fit_observable_curves;
  for (const SEXP &curve_sexp : fit_obj_rcpp) {
    Rcpp::List curve_rcpp(curve_sexp);
    const auto type = Rcpp::as<std::string>(curve_rcpp["type"]);
    if (type != "Gaussian")
      throw std::invalid_argument(
          "Each fit curve in fit.obj must have $type==\"Gaussian\"");
    const auto f = Rcpp::as<std::vector<Mean>>(curve_rcpp["f"]);
    const auto tau = Rcpp::as<SD>(curve_rcpp["tau"]);
    const Iota<Time> curve_times(f.size());
    auto mean_curve_gen = MakeFMapZipProductIterable(&std::make_pair<Time, const Mean &>, curve_times, f);
    MeanCurve mean_curve(mean_curve_gen.begin(), mean_curve_gen.end());
    ObservableCurve observable_curve = MakeObservableCurve(mean_curve, tau);
    fit_observable_curves.push_back(observable_curve);
  }

  // Rcpp::Rcout << "ZIP make result" << std::endl;
  auto &&xformed_curves_gen = ObservableCurveGenerator<ZipProducts>(
               fit_observable_curves,
               y_scale_baseline,
               curve_index_choices,
               peak_time_choices,
               x_shift_choices,
               x_scale_choices,
               sd_choices,
               sd_scale_choices,
               peak_height_choices,
               y_scale_choices
               );
  // Rcpp::Rcout << "make result gen initd" << std::endl;
  // const auto &log_weights_gen = LogWeightGenerator(
  // 						   observed_past_shrinkage_map,
  // 						   reasonable_future_shrinkage_map,
  // 						   n_future_neighbors,
  // 						   current_trajectory,
  // 						   historical_trajectories,
  // 						   xformed_curves_gen);
  // Rcpp::List curves_rcpp(xformed_curves_gen.size());
  // Rcpp::NumericVector log_weights_rcpp(xformed_curves_gen.size());
  Rcpp::NumericMatrix means(output_times.size(),xformed_curves_gen.size());
  Rcpp::NumericMatrix sds(output_times.size(),xformed_curves_gen.size());
  Rcpp::NumericVector log_weights(xformed_curves_gen.size());
  std::vector<Rcpp::DataFrame> curves_vec;
  std::vector<LogLikelihood> log_weights_vec;
  curves_vec.reserve(xformed_curves_gen.size());
  log_weights_vec.reserve(xformed_curves_gen.size());
  // std::ptrdiff_t curve_i = 0; // leads to ambiguity with []: ptr manip vs operator[]
  std::size_t curve_i = 0;
  for (const auto &xformed_curve : xformed_curves_gen) {
    // if (curve_i % 200 == 0) {
    //   Rcpp::Rcout << curve_i << std::endl;
    // }

    std::ptrdiff_t time_i = 0;
    for (Observable observable : InferencesForTimesGenerator(xformed_curve, output_times)) {
      means(time_i,curve_i) = observable.mean_;
      sds(time_i,curve_i) = observable.sd_;
      ++time_i;
    }
    log_weights[curve_i] = TrajectoryLogWeight(observed_past_shrinkage_map,
                                               reasonable_future_shrinkage_map,
                                               n_future_neighbors,
                                               current_trajectory,
                                               historical_trajectories,
                                               xformed_curve,
                                               bias_peaktime_mean,
                                               bias_peaktime_sd,
                                               bias_peaktime_shrinkage,
                                               bias_peakheight_mean,
                                               bias_peakheight_sd,
                                               bias_peakheight_shrinkage
                                               );
    ++curve_i;
  }
  Rcpp::List result_rcpp({{means, sds, log_weights}});
  return result_rcpp;
}

// todo bias shrinkage/scaling weight
// todo alternative approach: pseudocount prior + kernel conditional streaming
// todo try to use features of [[Rcpp::export]] to eliminate SEXP conversions when possible (typedefs complicate)

// todo clean and comment
// todo incremental, pseudocount prior + kernel-based conditional algorithm
// todo investigate huge compiled size

// clang-format off
// todo system includes
/* Local Variables: */
/* clang-format-style: "Google" */
/* flycheck-clang-language-standard: "c++14" */
/* flycheck-gcc-language-standard: "c++14" */
/* flycheck-clang-include-path: ("/usr/local/lib/R/site-library/Rcpp/include/" "/usr/share/R/include/") */
/* flycheck-gcc-include-path: ("/usr/local/lib/R/site-library/Rcpp/include/" "/usr/share/R/include/") */
/* End: */
// clang-format on
