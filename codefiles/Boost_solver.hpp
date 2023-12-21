//Find the super-adiabatic gradient that gives the given energy flux using Tom's Algorith from boost. Adapted from instructions available at:
// https://www.boost.org/doc/libs/1_74_0/libs/math/doc/html/math_toolkit/root_finding_examples/cbrt_eg.html

#include <boost/math/tools/roots.hpp>
#include "Fluxes.hpp"

template <class T> /// Functor for finding the super-adiabatic gradient from a flux
struct nabla_functor
{
    nabla_functor(T const & _L , T const & Re_crit , T const & A , T const & B , T const & l_nu , T const & K_cond , T const & nab_ad) : 
        L(_L) , F(Re_crit, A , B, l_nu , K_cond , nab_ad)
    {/* For fully molten/solid case, or no settling. Constructor stores values of coefficients and sets up Flux Structure*/ 
       
    }
    nabla_functor(T const & _L , T const & Re_crit , T const & A , T const & B , T const & l_nu , T const & K_cond , T const & nab_ad , T const& _DeltaS , T const &_phi, T const& Cp , T const & TdiffphiT , T const & PdiffphiP , T const & C  , char & _nabtype) : 
        L(_L) , F(Re_crit, A , B, l_nu , K_cond , nab_ad , Cp , TdiffphiT , PdiffphiP ) , DeltaS(_DeltaS) , phi(_phi) , Cgrav(C), nabtype(_nabtype)
    {/* For partial melt case. Constructor stores values of coefficients and sets up Flux Structure*/ }
    T operator()(T const& x) /// This operator is equal to zero if the super-adiabatic gradient gives the correct flux
    { 
        if(x < 1e-35){
          return 0;
        }
        T fx = F.LCond(x) + F.LConv(x) - L;
        if(nabtype == 'g'){
          fx = fx + F.LGrav(Cgrav , phi , DeltaS);
        }else
        if(nabtype == 'B'){
          fx = fx + F.LMix(x , DeltaS) + F.LGrav(Cgrav , phi , DeltaS);
        }
        //std::cout << F.FCond(x) << " " << Flux << std::endl;
        return fx;
    }
    private:
        T L  , DeltaS , phi , Cgrav;
        char nabtype;
        Fluxes<T> F; 
};

template <class T> /// Simpler case without gravitational settling (Curry et. al., 2023)
T boost_Tom(T L , T Re_crit , T A , T B , T l_nu , T K_cond , T nab_ad , T guess)
{
  // Solve using bracket_and_solve (no derivatives).
  using namespace boost::math::tools;           // For bracket_and_solve_root.

  T factor = 4;                                 // How big steps to take when searching.

  const boost::uintmax_t maxit = 50;            // Limit to maximum iterations.
  boost::uintmax_t it = maxit;                  // Initially our chosen max iterations, but updated with actual.
  bool is_rising = true;                        // So if result if f(guess) too low, then try increasing guess.
  int digits = std::numeric_limits<T>::digits;  // Maximum possible binary digits accuracy for type T.
  // Some fraction of digits is used to control how accurate to try to make the result.
  int get_digits = digits - 0;      
                         
  eps_tolerance<T> tol(get_digits);             // Set the tolerance.

  std::pair<T, T> r = bracket_and_solve_root(nabla_functor<T>(L, Re_crit , A , B , l_nu , K_cond , nab_ad), guess, factor, is_rising, tol, it);
  
  if (it >= maxit)
{
  std::cout << "Unable to locate solution with Brent's method in " << maxit << " iterations:"
    " Current best guess is between " << r.first << " and " << r.second << std::endl;
}

  return r.first + (r.second - r.first)/2;      // Midway between brackets is our result, if necessary we could
                                                // return the result as an interval here.                                               
}

template <class T> /// Case with gravitational settling flux (Abe)
T boost_Tom(T L , T Re_crit , T A , T B , T l_nu , T K_cond , T nab_ad , T DeltaS , T phi, T Cp , T TdiffphiT , T PdiffphiP , T C , T guess , char &nabtype)
{
  // Solve using bracket_and_solve (no derivatives).
  using namespace boost::math::tools;           // For bracket_and_solve_root.

  T factor = 4;                                 // How big steps to take when searching.

  const boost::uintmax_t maxit = 50;            // Limit to maximum iterations.
  boost::uintmax_t it = maxit;                  // Initially our chosen max iterations, but updated with actual.
  bool is_rising = true;                        // So if result if f(guess) too low, then try increasing guess.
  int digits = std::numeric_limits<T>::digits;  // Maximum possible binary digits accuracy for type T.
  // Some fraction of digits is used to control how accurate to try to make the result.
  int get_digits = digits - 0;      
                         
  eps_tolerance<T> tol(get_digits);             // Set the tolerance.
  //std::cout << "more complex solver" <<std::endl;
  std::pair<T, T> r = bracket_and_solve_root(nabla_functor<T>(L, Re_crit , A , B , l_nu , K_cond , nab_ad , DeltaS , phi , Cp , TdiffphiT , PdiffphiP , C , nabtype), guess, factor, is_rising, tol, it);
  //std::cout << "All values " << L << " " << Re_crit << " " << A << " " << B << " " << l_nu << " " << K_cond << " " << nab_ad << " " << DeltaS << " " << phi << " " << Cp << " " << TdiffphiT << " " << PdiffphiP << " " << C << " " << nabtype << std::endl;

  if (it >= maxit || (r.first < 1e-35 && r.first < 1e-35))
{
  std::cout << "Unable to locate solution with Brent's method in " << maxit << " iterations:"
    " Current best guess is between " << r.first << " and " << r.second << std::endl;
}
//std::cout << "sol " << r.first + (r.second - r.first)/2 << std::endl;
  return r.first + (r.second - r.first)/2;      // Midway between brackets is our result, if necessary we could
                                                // return the result as an interval here.                                               
}