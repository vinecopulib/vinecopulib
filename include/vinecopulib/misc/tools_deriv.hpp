// modified version of https://github.com/geotz/autodiff to make it C++11 
// compatible. Original license statement:

/*
MIT License

Copyright (c) 2017 George Tzoumas

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include <type_traits>
#include <cmath>

namespace std {
    template<class T>
    using decay_t = typename decay<T>::type;
    
    template<bool B, class T = void>
    using enable_if_t = typename enable_if<B, T>::type;
}

template<int N, typename... Ts> using NthTypeOf =
        typename std::tuple_element<N, std::tuple<Ts...>>::type;


// for debugging
template<class...> struct Foo;

template<class T> struct Expr {};

template<int K>
struct Const: Expr<Const<K>>
{
    using Value = std::integral_constant<int,K>;
    
    template<class...Var>
    static constexpr auto subs(Var...) -> int
    { 
        return Value::value; 
    }
    
    static constexpr auto diff() -> Const<K> 
    { 
        return Const{}; 
    }
    
    template<class V,class...Vs>
    static constexpr auto diff(V, Vs...) -> Const<0> 
    { 
        return Const<0>{}; 
    }
};

using Zero = Const<0>;
using Unit = Const<1>;

// template<int Acc>
// constexpr auto digit_parser() -> Const<Acc>
// {
//     return Const<Acc>{};
// }
// 
// // TODO: fix negative
// 
// template<int Acc,char D, char... Ds>
// constexpr auto digit_parser() -> Expr
// {
//     static_assert( D >= '0' && D <= '9', "invalid decimal digit" );
//     return digit_parser<Acc*10 + D-'0',Ds...>();
// }
// 
// template<char... Ds>
// constexpr auto operator ""_K() -> Expr
// {
//     return digit_parser<0,Ds...>();
// }

template<class T,int Id>
struct Lvalue
{
    T _val;
};

template<int Id>
struct Var: Expr<Var<Id>>
{
    using id_t = std::integral_constant<int,Id>;
    static constexpr auto id() -> id_t 
    { 
        return id_t{}; 
    }
    
    template<class T>
    constexpr auto operator()(T&& v) -> Lvalue<std::decay_t<T>,Id> const
    {
        return Lvalue<std::decay_t<T>,Id>{std::forward<T>(v)};
    }
    
    template<int...Is,class...Ts>
    static constexpr auto subs(const Lvalue<Ts,Is>&... vs)
        -> NthTypeOf<0, Ts...>
    {
        return Var::subs_multi( std::false_type{}, vs... );
    }
    
    static constexpr auto diff() -> Var<Id> { return Var{}; }
    
    template<class V,class... Vs>
    static constexpr auto diff(V, Vs... vs) 
    -> decltype(Const<V::id_t::value == Id>().diff(vs...))
    {
        using idv = typename V::id_t;
        return eval_diff(std::integral_constant<bool,idv::value == Id>{}).diff(vs...);
    }
    
private:
    static constexpr auto eval_diff( std::true_type ) -> Const<1>
    {
        return Const<1>{};
    }
    
    static constexpr auto eval_diff( std::false_type ) -> Const<0>
    {
        return Const<0>{};
    }
    
    // terminate only if var substituted
    static constexpr auto subs_multi( std::true_type ) -> decltype(0)
    {
        return 0;
    }
    
    template<bool B>
    static constexpr void subs_multi( std::integral_constant<bool,B> )
    {
        static_assert( B, "not all variables substituted" );
    }
    
    template<bool Found,int I1, int... Is,class T1,class... Ts>
    static constexpr auto subs_multi( 
        std::integral_constant<bool,Found>, 
        const Lvalue<T1,I1>& v1, 
        const Lvalue<Ts,Is>&... vs) 
    -> T1
    {
        return Var::subs_single<Found>(v1) + 
            Var::subs_multi(std::integral_constant<bool,Found || I1==Id>{}, vs...);
    }
    
    template<bool Found,class T>
    static constexpr auto subs_single( const Lvalue<T,Id>& v ) -> T
    {
        static_assert( !Found, "substitution already specified" );
        return v._val;
    }
    
    template<bool Found,int I,class T>
    static constexpr auto subs_single( const Lvalue<T,I>&  ) -> T
    {
        return T(0);
    }
};

template<class Op1,class Op2>
struct Add: Expr<Add<Op1,Op2>>
{
    template<int...Is,class...Ts>
    static constexpr auto subs(const Lvalue<Ts,Is>&... vs) 
        -> decltype(Op1::subs(vs...) + Op2::subs(vs...))
    {
        return Op1::subs(vs...) + Op2::subs(vs...);
    }

    static constexpr auto diff() -> Add { return Add{}; }

    template<class V,class... Vs>
    static constexpr auto diff(V v, Vs... vs) 
        -> decltype(
        (
            std::declval<decltype(Op1::diff(v))>() + std::declval<decltype(Op2::diff(v))>()
        ).diff(vs...)
    )
    {
        using D1 = decltype(Op1::diff(v));
        using D2 = decltype(Op2::diff(v));
        using Res = decltype((std::declval<D1>() + std::declval<D2>()).diff(vs...));
        return Res{};
    }
};

template<class Op>
struct Opp: Expr<Opp<Op>>
{
    template<int...Is,class...Ts>
    static constexpr auto subs(const Lvalue<Ts,Is>&... vs) 
        -> decltype(Op::subs(vs...))
    {
        return - Op::subs(vs...);
    }

    static constexpr auto diff() -> Opp<Op> { return Opp{}; }

    template<class V,class... Vs>
    static constexpr auto diff(V v, Vs... vs) 
        -> decltype((-std::declval<decltype(Op::diff(v))>()).diff(vs...))
    {
        using D = decltype(Op::diff(v));
        using Res = decltype((-std::declval<D>()).diff(vs...));
        return Res{};
    }
};

template<class Op1,class Op2>
struct Sub: Expr<Sub<Op1,Op2>>
{
    template<int...Is,class...Ts>
    static constexpr auto subs(const Lvalue<Ts,Is>&... vs) 
        -> decltype(Op1::subs(vs...) - Op2::subs(vs...)) 
    {
        return Op1::subs(vs...) - Op2::subs(vs...);
    }

    static constexpr auto diff() -> Sub<Op1, Op2> { return Sub{}; }

    template<class V,class... Vs>
    static constexpr auto diff(V v, Vs... vs)
    -> decltype((std::declval<decltype(Op1::diff(v))>() - std::declval<decltype(Op2::diff(v))>()).diff(vs...))
    {
        using D1 = decltype(Op1::diff(v));
        using D2 = decltype(Op2::diff(v));
        using Res = decltype((std::declval<D1>() - std::declval<D2>()).diff(vs...));
        return Res{};
    }
};

template<class Op1,class Op2>
struct Mul: Expr<Mul<Op1,Op2>>
{
    template<int...Is,class...Ts>
    static constexpr auto subs(const Lvalue<Ts,Is>&... vs)
        -> decltype(Op1::subs(vs...) * Op2::subs(vs...))
    {
        return Op1::subs(vs...) * Op2::subs(vs...);
    }

    static constexpr auto diff() -> Mul<Op1, Op2> { return Mul{}; }

    template<class V,class... Vs>
    static constexpr auto diff(V v, Vs... vs)
        -> decltype(
            (
                std::declval<decltype(Op1::diff(v))>() * std::declval<Op2>() +
                std::declval<Op1>() * std::declval<decltype(Op2::diff(v))>()
            ).diff(vs...)
        )
    {
        using D1 = decltype(Op1::diff(v));
        using D2 = decltype(Op2::diff(v));
        using Res = decltype((std::declval<D1>()*std::declval<Op2>() + std::declval<Op1>()*std::declval<D2>()).diff(vs...));
        std::cout << "called" << std::endl;
        return Res{};
    }
};

template<class Op,int K>
struct Pow: Expr<Pow<Op,K>>
{
    template<int...Is,class...Ts>
    static constexpr auto subs(const Lvalue<Ts,Is>&... vs)
        -> decltype(std::pow(Op::subs(vs...),K))
    {
        return std::pow(Op::subs(vs...),K);
    }

    static constexpr auto diff() -> Pow { return Pow{}; }

    template<class V,class... Vs>
    static constexpr auto diff(V v, Vs... vs)
    -> decltype(
        (
            Const<K>{}*(std::declval<Op>()^Const<K-1>{}) *
            std::declval<decltype(Op::diff(v))>()
        ).diff(vs...)
    )
    {
        using D = decltype(Op::diff(v));
        using Res = decltype((Const<K>{}*(std::declval<Op>()^Const<K-1>{})*std::declval<D>()).diff(vs...));
        return Res{};
    }
};

template<class>
struct IsPow: std::false_type
{};

template<class T,int K>
struct IsPow<Expr<Pow<T,K>>>: std::true_type
{};

template<class Op1,class Op2>
struct Div: Expr<Div<Op1,Op2>>
{
    template<int...Is,class...Ts>
    static constexpr auto subs(const Lvalue<Ts,Is>&... vs)
    -> decltype(Op1::subs(vs...) / Op2::subs(vs...))
    {
        return Op1::subs(vs...) / Op2::subs(vs...);
    }

    static constexpr auto diff() -> Div { return Div{}; }

    template<class V,class... Vs>
    static constexpr auto diff(V v, Vs... vs)
    -> decltype(
        (
            (
                std::declval<decltype(Op1::diff(v))>() * std::declval<Op2>() -
                std::declval<Op1>() * std::declval<decltype(Op2::diff(v))>()
            ) /
            (std::declval<Op2>()^std::declval<Const<2>>())
        ).diff(vs...)
        )
    {
        using D1 = decltype(Op1::diff(v));
        using D2 = decltype(Op2::diff(v));
        using Res = decltype(((std::declval<D1>()*std::declval<Op2>() -
                std::declval<Op1>()*std::declval<D2>())/
                    (std::declval<Op2>()^std::declval<Const<2>>())).diff(vs...));
        return Res{};
    }
};

//template<class> struct Cos;

template<class Op>
struct Sin: Expr<Sin<Op>>
{
    template<int...Is,class...Ts>
    static constexpr auto subs(const Lvalue<Ts,Is>&... vs)
    -> decltype(std::sin(Op::subs(vs...)))
    {
        return std::sin(Op::subs(vs...));
    }

    static constexpr auto diff() -> Sin { return Sin{}; }

    template<class V,class... Vs>
    static constexpr auto diff(V v, Vs... vs)
    -> decltype(
        (cos(std::declval<Op>())*std::declval<decltype(Op::diff(v))>()).diff(vs...)
    )
    {
        using D = decltype(Op::diff(v));
        using Res = decltype((cos(std::declval<Op>())*std::declval<D>()).diff(vs...));
        return Res{};
    }
};

template<class Op>
struct Cos: Expr<Cos<Op>>
{
    template<int...Is,class...Ts>
    static constexpr auto subs(const Lvalue<Ts,Is>&... vs)
    -> decltype(std::cos(Op::subs(vs...)))
    {
        return std::cos(Op::subs(vs...));
    }

    static constexpr auto diff() -> Cos { return Cos{}; }

    template<class V,class... Vs>
    static constexpr auto diff(V v, Vs... vs)
    -> decltype(
        (
            Const<-1>{}*sin(std::declval<Op>())*std::declval<decltype(Op::diff(v))>()
        ).diff(vs...)
    )
    {
        using D = decltype(Op::diff(v));
        using Res = decltype((Const<-1>{}*sin(std::declval<Op>())*std::declval<D>()).diff(vs...));
        return Res{};
    }
};

// --------------------------------
//   BASIC operators
// --------------------------------
template<class Op1,class Op2>
constexpr auto operator+(Expr<Op1>, Expr<Op2>) -> Add<Op1, Op2>
{ 
    return Add<Op1,Op2>{}; 
}

template<class Op1,class Op2>
constexpr auto operator-(Expr<Op1>, Expr<Op2>) 
-> decltype(Op1{} + Const<-1>{}*Op2{})
{ 
    return Op1{} + Const<-1>{}*Op2{}; 
}

template<class Op>
constexpr auto operator-(Expr<Op>) -> decltype(Const<-1>{}*Op{})
{ 
    return Const<-1>{}*Op{}; 
}

template<class Op1,class Op2>
constexpr auto operator*(Expr<Op1>, Expr<Op2>) -> Mul<Op1, Op2>
{ 
    return Mul<Op1,Op2>{}; 
}

template<class Op,int K>
constexpr auto operator^(Expr<Op>, Expr<Const<K>>) -> Pow<Op,K>
{ 
    return Pow<Op,K>{}; 
}

template<class Op1,class Op2>
constexpr auto operator/(Expr<Op1>, Expr<Op2>) -> Div<Op1,Op2>
{ 
    return Div<Op1,Op2>{}; 
}

template<class Op>
constexpr auto sin(Expr<Op>) -> Sin<Op>
{ 
    return Sin<Op>{}; 
}

template<class Op>
constexpr auto cos(Expr<Op>) -> Cos<Op>
{ 
    return Cos<Op>{}; 
}

// --------------------------------
//   ADD trivial simplifications
// --------------------------------

template<class Op>
constexpr auto operator+(Expr<Op>, Expr<Op>) -> decltype(Const<2>{}*Op{})
{ 
    return Const<2>{}*Op{}; 
}

template<int K,class Op>
constexpr auto operator+(Expr<Mul<Const<K>,Op>>, Expr<Op>) 
    -> decltype(Const<K+1>{}*Op{})
{ 
    return Const<K+1>{}*Op{}; 
}

template<int K,class Op>
constexpr auto operator+(Expr<Op>, Expr<Mul<Const<K>,Op>>)
    -> decltype(Const<K+1>{}*Op{})
{ 
    return Const<K+1>{}*Op{}; 
}

template<class Op>
constexpr auto operator+(Expr<Zero>, Expr<Op>) -> Op
{ 
    return Op{}; 
}

template<class Op>
constexpr auto operator+(Expr<Op>,Expr<Zero>) -> Op
{ 
    return Op{}; 
}

constexpr auto operator+(Expr<Zero>,Expr<Zero>) -> Zero
{ 
    return Zero{}; 
}

template<int K,class=std::enable_if_t<K!=0>>
constexpr auto operator+(Expr<Zero>,Expr<Const<K>>) -> Const<K>
{ 
    return Const<K>{}; 
}

template<int K,class=std::enable_if_t<K!=0>>
constexpr auto operator+(Expr<Const<K>>, Expr<Zero>) -> Const<K>
{ 
    return Const<K>{}; 
}

template<int K1,int K2,class=std::enable_if_t<K1!=0 && K2!=0 && K1!=K2>>
constexpr auto operator+(Expr<Const<K1>>, Expr<Const<K2>>) -> Const<K1+K2>
{ 
    return Const<K1+K2>{}; 
}

//// --------------------------------
////   OPP (unary minus) trivial simplifications
//// --------------------------------

template<int K>
constexpr auto operator-(Expr<Const<K>>) -> Const<-K>
{ 
    return Const<-K>{}; 
}

// --------------------------------
//   SUB trivial simplifications
// --------------------------------

template<class Op>
constexpr auto operator-(Expr<Op>,Expr<Op>) -> Zero
{ 
    return Zero{}; 
}

template<class Op>
constexpr auto operator-(Expr<Zero>, Expr<Op>) -> decltype(Const<-1>{} * Op{})
{ 
    return Const<-1>{} * Op{}; 
}

template<class Op>
constexpr auto operator-(Expr<Op>,Expr<Zero>) -> Op
{ 
    return Op{}; 
}

constexpr auto operator-(Expr<Zero>,Expr<Zero>) -> Zero
{ 
    return Zero{}; 
}

//template<int K,class=std::enable_if_t<K!=0>>
//constexpr auto operator-(Expr<Zero>,Expr<Const<K>>) { return Const<-K>{}; }

//template<int K,class=std::enable_if_t<K!=0>>
//constexpr auto operator-(Expr<Const<K>>, Expr<Zero>) { return Const<K>{}; }

//template<int K1,int K2,class=std::enable_if_t<K1!=0 && K2!=0>>
//constexpr auto operator-(Expr<Const<K1>>, Expr<Const<K2>>) { return Const<K1-K2>{}; }

// --------------------------------
//   MUL trivial simplifications
// --------------------------------

template<class Op,class=std::enable_if_t<!IsPow<Op>::value>>
constexpr auto operator*(Expr<Op>, Expr<Op>) -> decltype(Op{}^Const<2>{})
{ 
    return Op{}^Const<2>{}; 
}

template<int K1,int K2,class Op,class=std::enable_if_t<IsPow<Op>::value>>
constexpr auto operator*(Expr<Pow<Op,K1>>, Expr<Pow<Op,K2>>) 
    -> decltype(Op{}^Const<K1+K2>{})
{ 
    return Op{}^Const<K1+K2>{}; 
}

template<int K1,int K2,class Op>
constexpr auto operator*(Const<K1>, Expr<Mul<Const<K2>,Op>>) 
    -> decltype(Const<K1*K2>{}*Op{})
{ 
    return Const<K1*K2>{}*Op{}; 
}

template<int K1,int K2,class Op>
constexpr auto operator*(Expr<Mul<Const<K2>,Op>>, Const<K1>) 
    -> decltype(Const<K1*K2>{}*Op{})
{ 
    return Const<K1*K2>{}*Op{}; 
}

template<class Op>
constexpr auto operator*(Expr<Zero>, Expr<Op>) -> Zero
{ 
    return Zero{}; 
}

template<class Op>
constexpr auto operator*(Expr<Op>,Expr<Zero>) 
    -> Zero
{ 
    return Zero{}; 
}

template<class Op>
constexpr auto operator*(Expr<Unit>, Expr<Op>) -> Op
{ 
    return Op{}; 
}

template<class Op>
constexpr auto operator*(Expr<Op>,Expr<Unit>) -> Op 
{ 
    return Op{}; 
}

constexpr auto operator*(Expr<Unit>, Expr<Zero>) -> Zero
{ 
    return Zero{}; 
}

constexpr auto operator*(Expr<Zero>, Expr<Unit>) -> Zero
{ 
    return Zero{}; 
}

constexpr auto operator*(Expr<Zero>, Expr<Zero>) -> Zero
{ 
    return Zero{}; 
}

constexpr auto operator*(Expr<Unit>, Expr<Unit>) -> Unit
{ 
    return Unit{}; 
}

template<int K1,int K2,class=std::enable_if_t<K1!=0 && K1!=1 && K2!=0 && K2!=1>>
constexpr auto operator*(Expr<Const<K1>>, Expr<Const<K2>>) -> Const<K1*K2>
{ 
    return Const<K1*K2>{}; 
}

// --------------------------------
//   POW trivial simplifications
// --------------------------------

template<class Op>
constexpr auto operator^(Expr<Op>, Expr<Unit>) -> Op
{ 
    return Op{}; 
}

template<class Op>
constexpr auto operator^(Expr<Op>, Expr<Zero>)  -> Unit
{ 
    return Unit{}; 
}

template<int K,int M,class Op,class=std::enable_if_t<(M>1)>>
constexpr auto operator^(Expr<Pow<Op,K>>, Expr<Const<M>>) -> Pow<Op, K*M>
{ 
    return Pow<Op,K*M>{}; 
}

// --------------------------------
//   DIV trivial simplifications
// --------------------------------

template<class Op>
constexpr auto operator/(Expr<Zero>, Expr<Op>) -> Zero
{ 
    return Zero{}; 
}

template<class Op,int K=0>
constexpr void operator/(Expr<Op>,Expr<Zero>) 
{ 
    static_assert( K!=0, "division by zero" ); 
}

template<int K=0>
constexpr void operator/(Expr<Zero>,Expr<Zero>) 
{ 
    static_assert( K!=0, "NaN (0/0)" ); 
}

template<class Op>
constexpr auto operator/(Expr<Op>,Expr<Unit>) -> Op
{ 
    return Op{}; 
}

constexpr auto operator/(Expr<Zero>, Expr<Unit>) -> Zero
{ 
    return Zero{}; 
}
constexpr auto operator/(Expr<Unit>, Expr<Unit>) -> Unit
{ 
    return Unit{}; 
}

template<class Op>
constexpr auto operator/(Expr<Op>,Expr<Op>) -> Unit
{ 
    return Unit{}; 
}

template<int K1,int K2,int GCD>
constexpr auto simplify_frac(Expr<Const<K1>>, Expr<Const<K2>>, std::integral_constant<int,GCD>)
-> Div<Const<K1/GCD>,Const<K2/GCD>>
{
    return Div<Const<K1/GCD>,Const<K2/GCD>>{};
}

template<int K1,int K2>
constexpr auto simplify_frac(Expr<Const<K1>>, Expr<Const<K2>>, std::integral_constant<int,K2>)
-> Const<K1/K2>
{
    return Const<K1/K2>{};
}

template<int K1,int K2,class=std::enable_if_t<K1!=0 && K2!=0 && K2!=1 && K1!=K2>>
constexpr auto operator/(Expr<Const<K1>> a, Expr<Const<K2>> b) 
-> decltype(simplify_frac(a, b, gcd( std::integral_constant<int,K1>{}, std::integral_constant<int,K2>{} )))
{
    auto g = gcd( std::integral_constant<int,K1>{}, std::integral_constant<int,K2>{} );
    return simplify_frac(a, b, g);
}

// int main()
// {
//     using namespace std;
//     using X = Var<0>;
//     using Y = Var<1>;
//     using Z = Var<2>;
//     X x;
//     Y y;
//     Z z;
//     auto r1 = (x*y) / (x+y);
//     auto r2 = sin(x)*sin(x) + cos(x)*cos(x);
//     auto r3 = (Const<3>{} / Const<12>{});// (x+y)*(x+y);
//  
//     // this works
//     cerr << r1.subs( x(2.), y(3.)) << endl;
//     cerr << r1.diff(x).subs( x(2.), y(3.)) << endl;
//     cerr << r1.diff(x,y).subs( x(2.), y(3.)) << endl;
//      
//     // this too
//     cerr << r2.subs( x(2.), y(3.)) << endl;
//     cerr << r2.diff(x).subs( x(2.), y(3.)) << endl;
//     cerr << r2.diff(x,y).subs( x(2.), y(3.)) << endl;
// 
//     // this doesn't
//     cerr << r3.subs( x(2.), y(3.)) << endl;
//     cerr << r3.diff(x).subs( x(2.), y(3.)) << endl;
//     cerr << r3.diff(x,y).subs( x(2.), y(3.)) << endl;
// }
