#ifndef ACE_RADIAL_FUNCTIONS_H
#define ACE_RADIAL_FUNCTIONS_H

#include "ace_arraynd.h"
#include "ace_types.h"
#include <functional>

using namespace std;

//typedef void (*RadialFunctions)(DOUBLE_TYPE x);
typedef std::function<void(DOUBLE_TYPE)> RadialFunctions;

/**
 * Class that implement spline interpolation and caching for radial functions
 */
class SplineInterpolator {
public:
    DOUBLE_TYPE cutoff = 0; ///< cutoff
    int ntot = 10000; ///< Number of bins for look-up tables.
    int nlut; ///< number of nodes in look-up table
    DOUBLE_TYPE invrscalelookup; ///< inverse of conversion coefficient from distance to lookup table within cutoff range
    DOUBLE_TYPE rscalelookup; ///< conversion coefficient from distance to lookup table within cutoff range
    int num_of_functions;///< number of functions to spline-interpolation

    Array1D<DOUBLE_TYPE> values = Array1D<DOUBLE_TYPE>("values"); ///< shape: [func_ind]
    Array1D<DOUBLE_TYPE> derivatives = Array1D<DOUBLE_TYPE>("derivatives");///< shape: [func_ind]
    Array3D<DOUBLE_TYPE> lookupTable = Array3D<DOUBLE_TYPE>("lookupTable");///< shape: [ntot+1][func_ind][4]

    /**
     * Setup splines
     *
     * @param num_of_functions number of functions
     * @param func subroutine, that update `values` and `dvalues` arrays
     * @param values values
     * @param dvalues derivatives
     */
    void setupSplines(int num_of_functions, RadialFunctions func,
                      DOUBLE_TYPE *values,
                      DOUBLE_TYPE *dvalues, int ntot, DOUBLE_TYPE cutoff);

    /**
     * Populate `values` and `derivatives` arrays with a spline-interpolation for
     * all functions
     *
     * @param r
     *
     * @return: populate 'values' and 'derivatives'
     */
    void calcSplines(DOUBLE_TYPE r);

    /**
     * Populate `values` and `derivatives` arrays with a spline-interpolation for
     * all functions
     *
     * @param r
     *
     * @return: populate 'values' and 'derivatives'
     */
    void calcSplines(DOUBLE_TYPE r, SHORT_INT_TYPE func_ind);
};


class AbstractRadialBasis {
public:
    SPECIES_TYPE nelements = 0; ///< number of elements
    Array2D<DOUBLE_TYPE> cut = Array2D<DOUBLE_TYPE>("cut"); ///< cutoffs, shape: [nelements][nelements]
    Array2D<DOUBLE_TYPE> dcut = Array2D<DOUBLE_TYPE>("dcut"); ///< decay of cutoff, shape: [nelements][nelements]
    DOUBLE_TYPE cutoff = 0; ///< cutoff

    int ntot = 10000; ///< Number of bins for look-up tables.
    LS_TYPE lmax = 0; ///< maximum value of `l`
    NS_TYPE nradial = 0;  ///< maximum number `n` of radial functions \f$ R_{nl}(r) \f$
    NS_TYPE nradbase = 0; ///< number of radial basis functions \f$ g_k(r) \f$

    // Arrays for look-up tables.
    Array2D<SplineInterpolator> splines_gk; ///< array of spline interpolator to store g_k, shape: [nelements][nelements]
    Array2D<SplineInterpolator> splines_rnl; ///< array of spline interpolator to store R_nl, shape: [nelements][nelements]
    Array2D<SplineInterpolator> splines_hc; ///< array of spline interpolator to store R_nl shape: [nelements][nelements]
    //--------------------------------------------------------------------------

    string radbasename = "ChebExpCos"; ///< type of radial basis functions \f$ g_{k}(r) \f$ (default="ChebExpCos")

    /**
   Arrays to store radial functions.
   */
    Array1D<DOUBLE_TYPE> gr = Array1D<DOUBLE_TYPE>("gr"); ///< g_k(r) functions, shape: [nradbase]
    Array1D<DOUBLE_TYPE> dgr = Array1D<DOUBLE_TYPE>("dgr"); ///< derivatives of g_k(r) functions, shape: [nradbase]
    Array2D<DOUBLE_TYPE> fr = Array2D<DOUBLE_TYPE>("fr");  ///< R_nl(r) functions, shape: [nradial][lmax+1]
    Array2D<DOUBLE_TYPE> dfr = Array2D<DOUBLE_TYPE>(
            "dfr"); ///< derivatives of R_nl(r) functions, shape: [nradial][lmax+1]



    DOUBLE_TYPE cr; ///< hard-core repulsion
    DOUBLE_TYPE dcr; ///< derivative of hard-core repulsion

    Array5D<DOUBLE_TYPE> crad = Array5D<DOUBLE_TYPE>(
            "crad"); ///< expansion coefficients of radial functions into radial basis function, see Eq. (27) of PRB, shape:  [nelements][nelements][lmax + 1][nradial][nradbase]
    Array2D<DOUBLE_TYPE> lambda = Array2D<DOUBLE_TYPE>(
            "lambda"); ///< distance scaling parameter Eq.(24) of PRB,  shape: [nelements][nelements]

    Array2D<DOUBLE_TYPE> prehc = Array2D<DOUBLE_TYPE>(
            "prehc"); ///< hard-core repulsion coefficients (prefactor), shape: [nelements][nelements]
    Array2D<DOUBLE_TYPE> lambdahc = Array2D<DOUBLE_TYPE>(
            "lambdahc");; ///< hard-core repulsion coefficients (lambdahc), shape: [nelements][nelements]

    virtual void
    evaluate(DOUBLE_TYPE r, NS_TYPE nradbase_c, NS_TYPE nradial_c, SPECIES_TYPE mu_i, SPECIES_TYPE mu_j) = 0;

    virtual void
    init(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, int ntot, SPECIES_TYPE nelements, DOUBLE_TYPE cutoff,
         string radbasename = "ChebExpCos") = 0;

    /**
    * Function that sets up the look-up tables for spline-representation of radial functions.
    */
    virtual void setuplookupRadspline() = 0;

    virtual AbstractRadialBasis *clone() const = 0;

    virtual ~AbstractRadialBasis() = default;
};


/**
Class to store radial functions and their associated functions. \n
*/
class ACERadialFunctions final : public AbstractRadialBasis {
public:

    //--------------------------------------------------------------------------

    /**
    Arrays to store Chebyshev polynomials.
    */
    Array1D<DOUBLE_TYPE> cheb = Array1D<DOUBLE_TYPE>(
            "cheb"); ///< Chebyshev polynomials of the first kind, shape: [nradbase+1]
    Array1D<DOUBLE_TYPE> dcheb = Array1D<DOUBLE_TYPE>(
            "dcheb"); ///< derivatives Chebyshev polynomials of the first kind, shape: [nradbase+1]
    Array1D<DOUBLE_TYPE> cheb2 = Array1D<DOUBLE_TYPE>(
            "cheb2"); ///< Chebyshev polynomials of the second kind, shape: [nradbase+1]

    //--------------------------------------------------------------------------

    /**
     * Default constructor
     */
    ACERadialFunctions() = default;

    /**
     * Parametrized constructor
     *
     * @param nradb number of radial basis function \f$ g_k(r) \f$ - nradbase
     * @param lmax maximum orbital moment - lmax
     * @param nradial maximum n-index of radial functions  \f$ R_{nl}(r) \f$ - nradial
     * @param ntot   Number of bins for spline look-up tables.
     * @param nelements   numer of elements
     * @param cutoff cutoff
     * @param radbasename  type of radial basis function \f$ g_k(r) \f$ (default: "ChebExpCos")
     */
    ACERadialFunctions(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, int ntot, SPECIES_TYPE nelements,
                       DOUBLE_TYPE cutoff, string radbasename = "ChebExpCos");

    /**
     * Initialize arrays for given parameters
     *
     * @param nradb number of radial basis function \f$ g_k(r) \f$ - nradbase
     * @param lmax maximum orbital moment - lmax
     * @param nradial maximum n-index of radial functions  \f$ R_{nl}(r) \f$ - nradial
     * @param ntot   Number of bins for spline look-up tables.
     * @param nelements   numer of elements
     * @param cutoff cutoff
     * @param radbasename  type of radial basis function \f$ g_k(r) \f$ (default: "ChebExpCos")
     */
    void init(NS_TYPE nradb, LS_TYPE lmax, NS_TYPE nradial, int ntot, SPECIES_TYPE nelements, DOUBLE_TYPE cutoff,
              string radbasename = "ChebExpCos") override;

    /**
     * Destructor
     */
    ~ACERadialFunctions() override = default;

    /**
    * Function that computes Chebyshev polynomials of first and second kind
    * to setup the radial functions and the derivatives
    *
    * @param n  maximum polynom order
    * @param x
    *
    * @returns fills cheb, dcheb and cheb2 arrays
    */
    void calcCheb(NS_TYPE n, DOUBLE_TYPE x);

    /**
     * Function that computes radial basis functions  \$f g_k(r) \$f, see Eq.(21) of PRB paper
     * @param lam \$f \lambda \$f  parameter, see eq. (24) of PRB paper
     * @param cut cutoff
     * @param dcut cutoff decay
     * @param r distance
     *
     * @return  function fills gr and dgr arrays
     */
    void radbase(DOUBLE_TYPE lam, DOUBLE_TYPE cut, DOUBLE_TYPE dcut, DOUBLE_TYPE r);

    /**
     *   Function that computes radial core repulsion \$f f_{core} = pre \exp( - \lambda r^2 ) / r \$f,
     *   and its derivative, see Eq.(27) of implementation notes.
     *
     * @param r  distance
     * @param pre prefactor: read from input, depends on pair of atoms mu_i mu_j
     * @param lambda exponent: read from input, depends on pair of atoms mu_i mu_j
     * @param cutoff cutoff distance: read from input, depends on pair of atoms mu_i mu_j
     * @param cr  (out) hard core repulsion
     * @param dcr (out) derivative of hard core repulsion
     */
    static void radcore(DOUBLE_TYPE r, DOUBLE_TYPE pre, DOUBLE_TYPE lambda, DOUBLE_TYPE cutoff, DOUBLE_TYPE &cr,
                        DOUBLE_TYPE &dcr);

    /**
     * Function that sets up the look-up tables for spline-representation of radial functions.
     */
    void setuplookupRadspline() override;

    /**
     * Function that computes radial functions \f$ R_{nl}(r)\f$  (see Eq. 27 from PRB paper)
     * and its derivatives for all range of n,l,
     * ONLY if radial basis functions (gr and dgr) are computed.
     * @param elei first species type
     * @param elej second species type
     *
     * @return fills in fr, dfr arrays
     */
    void radfunc(SPECIES_TYPE elei, SPECIES_TYPE elej);

    /**
     * Compute all radial functions R_nl(r), radial basis functions g_k(r) and hard-core repulsion function hc(r)
     *
     * @param r distance
     * @param nradbase_c
     * @param nradial_c
     * @param mu_i
     * @param mu_j
     *
     * @return update gr(k), dgr(k), fr(n,l), dfr(n,l), cr, dcr
     */
    void evaluate(DOUBLE_TYPE r, NS_TYPE nradbase_c, NS_TYPE nradial_c, SPECIES_TYPE mu_i, SPECIES_TYPE mu_j) override;

    void chebExpCos(DOUBLE_TYPE lam, DOUBLE_TYPE cut, DOUBLE_TYPE dcut, DOUBLE_TYPE r);

    void chebPow(DOUBLE_TYPE lam, DOUBLE_TYPE cut, DOUBLE_TYPE dcut, DOUBLE_TYPE r);

    void chebLinear(DOUBLE_TYPE lam, DOUBLE_TYPE cut, DOUBLE_TYPE dcut, DOUBLE_TYPE r);

    /**
     * Setup all radial functions for element pair mu_i-mu_j and distance r
     * @param mu_i first specie type
     * @param mu_j second specie type
     * @param r distance
     * @return update fr(nr, l),  dfr(nr, l)
     */
    void all_radfunc(SPECIES_TYPE mu_i, SPECIES_TYPE mu_j, DOUBLE_TYPE r);

    ACERadialFunctions *clone() const override {
        return new ACERadialFunctions(*this);
    };
};

#endif