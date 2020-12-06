#ifndef LaserTherm_Utils_Mixins_hpp
#define LaserTherm_Utils_Mixins_hpp

/** @file Mixins.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/02/19
 */

template<typename REAL>
class MixinBase
{
 public:
  using real_type = REAL;
};

#define MAKE_ADD_MEMBER_MIXIN(NAME)                                        \
  template<typename BASE>                                                  \
  class Add##NAME : public BASE                                            \
  {                                                                        \
   public:                                                                 \
    using real_type = typename BASE::real_type;                            \
                                                                           \
   protected:                                                              \
    real_type NAME;                                                        \
                                                                           \
   public:                                                                 \
    void      set##NAME(const real_type& a) { NAME = a; }                  \
    real_type get##NAME() const { return NAME; }                           \
  };                                                                       \
  template<typename T>                                                     \
  constexpr auto has##NAME(T x)->decltype(x.get##NAME(), std::true_type{}) \
  {                                                                        \
    return {};                                                             \
  }                                                                        \
  constexpr auto has##NAME(...)->std::false_type { return {}; }

#define MAKE_ADD_OPTIONAL_MEMBER_MIXIN(NAME)                             \
  template<typename BASE>                                                \
  class AddOptional##NAME : public BASE                                  \
  {                                                                      \
   public:                                                               \
    using real_type = typename BASE::real_type;                          \
                                                                         \
   protected:                                                            \
    std::optional<real_type> NAME;                                       \
                                                                         \
   public:                                                               \
    void                     set##NAME(const real_type& a) { NAME = a; } \
    std::optional<real_type> get##NAME() const { return NAME; }          \
  };                                                                     \
  template<typename T>                                                   \
  constexpr auto hasOptional##NAME(T x)->decltype(x.get##NAME(),         \
                                                  std::true_type{})      \
  {                                                                      \
    return {};                                                           \
  }                                                                      \
  constexpr auto hasOptional##NAME(...)->std::false_type { return {}; }

#endif  // include protector
