
// Code generated by stanc f556d0d
#include <stan/model/model_header.hpp>
namespace dgp_lotka_volterra_model_namespace {

template <typename T, typename S>
std::vector<T> resize_to_match__(std::vector<T>& dst, const std::vector<S>& src) {
  dst.resize(src.size());
  return dst;
}

template <typename T>
Eigen::Matrix<T, -1, -1>
resize_to_match__(Eigen::Matrix<T, -1, -1>& dst, const Eigen::Matrix<T, -1, -1>& src) {
  dst.resize(src.rows(), src.cols());
  return dst;
}

template <typename T>
Eigen::Matrix<T, 1, -1>
resize_to_match__(Eigen::Matrix<T, 1, -1>& dst, const Eigen::Matrix<T, 1, -1>& src) {
  dst.resize(src.size());
  return dst;
}

template <typename T>
Eigen::Matrix<T, -1, 1>
resize_to_match__(Eigen::Matrix<T, -1, 1>& dst, const Eigen::Matrix<T, -1, 1>& src) {
  dst.resize(src.size());
  return dst;
}
std::vector<double> to_doubles__(std::initializer_list<double> x) {
  return x;
}

std::vector<stan::math::var> to_vars__(std::initializer_list<stan::math::var> x) {
  return x;
}

inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}

inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}


using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math; 

static int current_statement__ = 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 29, column 2 to column 27)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 30, column 2 to column 28)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 31, column 2 to column 27)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 32, column 2 to line 35, column 42)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 38, column 2 to column 17)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 39, column 2 to column 26)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 41, column 2 to column 32)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 42, column 2 to column 32)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 43, column 2 to column 36)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 44, column 2 to column 36)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 45, column 2 to column 34)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 46, column 2 to column 34)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 47, column 2 to column 40)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 48, column 2 to column 40)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 51, column 4 to column 56)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 52, column 4 to column 51)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 50, column 17 to line 54, column 3)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 50, column 2 to line 54, column 3)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 22, column 2 to column 19)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 23, column 2 to column 13)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 7, column 4 to column 18)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 8, column 4 to column 18)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 10, column 4 to column 26)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 11, column 4 to column 25)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 12, column 4 to column 26)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 13, column 4 to column 26)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 15, column 4 to column 40)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 16, column 4 to column 42)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 17, column 4 to column 28)",
                                                      " (in '/Users/hyunjimoon/Dropbox/stan/casestudy/lotka-volterra/dgp_lotka-volterra.stan', line 6, column 26 to line 18, column 3)"};


template <typename T0__, typename T1__, typename T2__, typename T3__>
std::vector<typename boost::math::tools::promote_args<T0__, T1__, T2__,
T3__>::type>
dz_dt(const T0__& t, const std::vector<T1__>& z,
      const std::vector<T2__>& theta, const std::vector<T3__>& x_r,
      const std::vector<int>& x_i, std::ostream* pstream__) {
  using local_scalar_t__ = typename boost::math::tools::promote_args<T0__,
          T1__,
          T2__,
          T3__>::type;
  const static bool propto__ = true;
  (void) propto__;
  
  try {
    local_scalar_t__ u;
    
    current_statement__ = 21;
    u = std::numeric_limits<double>::quiet_NaN();
    current_statement__ = 21;
    u = z[(1 - 1)];
    local_scalar_t__ v;
    
    current_statement__ = 22;
    v = std::numeric_limits<double>::quiet_NaN();
    current_statement__ = 22;
    v = z[(2 - 1)];
    local_scalar_t__ alpha;
    
    current_statement__ = 23;
    alpha = std::numeric_limits<double>::quiet_NaN();
    current_statement__ = 23;
    alpha = theta[(1 - 1)];
    local_scalar_t__ beta;
    
    current_statement__ = 24;
    beta = std::numeric_limits<double>::quiet_NaN();
    current_statement__ = 24;
    beta = theta[(2 - 1)];
    local_scalar_t__ gamma;
    
    current_statement__ = 25;
    gamma = std::numeric_limits<double>::quiet_NaN();
    current_statement__ = 25;
    gamma = theta[(3 - 1)];
    local_scalar_t__ delta;
    
    current_statement__ = 26;
    delta = std::numeric_limits<double>::quiet_NaN();
    current_statement__ = 26;
    delta = theta[(4 - 1)];
    local_scalar_t__ du_dt;
    
    current_statement__ = 27;
    du_dt = std::numeric_limits<double>::quiet_NaN();
    current_statement__ = 27;
    du_dt = ((alpha - (beta * v)) * u);
    local_scalar_t__ dv_dt;
    
    current_statement__ = 28;
    dv_dt = std::numeric_limits<double>::quiet_NaN();
    current_statement__ = 28;
    dv_dt = ((-gamma + (delta * u)) * v);
    current_statement__ = 29;
    return stan::math::array_builder<local_scalar_t__>().add(du_dt)
        .add(dv_dt).array();
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
  }
  
}

struct dz_dt_functor__ {
template <typename T0__, typename T1__, typename T2__, typename T3__>
std::vector<typename boost::math::tools::promote_args<T0__, T1__, T2__,
T3__>::type>
operator()(const T0__& t, const std::vector<T1__>& z,
           const std::vector<T2__>& theta, const std::vector<T3__>& x_r,
           const std::vector<int>& x_i, std::ostream* pstream__)  const 
{
return dz_dt(t, z, theta, x_r, x_i, pstream__);
}
};

class dgp_lotka_volterra_model : public model_base_crtp<dgp_lotka_volterra_model> {

 private:
  int pos__;
  int N;
  std::vector<double> ts;
 
 public:
  ~dgp_lotka_volterra_model() { }
  
  std::string model_name() const { return "dgp_lotka_volterra_model"; }
  
  dgp_lotka_volterra_model(stan::io::var_context& context__,
                           unsigned int random_seed__ = 0,
                           std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    typedef double local_scalar_t__;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "dgp_lotka_volterra_model_namespace::dgp_lotka_volterra_model";
    (void) function__;  // suppress unused var warning
    
    try {
      
      pos__ = 1;
      context__.validate_dims("data initialization","N","int",
          context__.to_vec());
      
      current_statement__ = 19;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 20;
      validate_non_negative_index("ts", "N", N);
      context__.validate_dims("data initialization","ts","double",
          context__.to_vec(N));
      ts = std::vector<double>(N, 0);
      
      current_statement__ = 20;
      assign(ts, nil_index_list(), context__.vals_r("ts"),
        "assigning variable ts");
      current_statement__ = 19;
      current_statement__ = 19;
      check_greater_or_equal(function__, "N", N, 0);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename T__>
  T__ log_prob(std::vector<T__>& params_r__, std::vector<int>& params_i__,
               std::ostream* pstream__ = 0) const {
    typedef T__ local_scalar_t__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "dgp_lotka_volterra_model_namespace::log_prob";
(void) function__;  // suppress unused var warning

    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    
    try {
      
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob() 
    
  template <typename RNG>
  void write_array(RNG& base_rng__, std::vector<double>& params_r__,
                   std::vector<int>& params_i__, std::vector<double>& vars__,
                   bool emit_transformed_parameters__ = true,
                   bool emit_generated_quantities__ = true,
                   std::ostream* pstream__ = 0) const {
    typedef double local_scalar_t__;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "dgp_lotka_volterra_model_namespace::write_array";
(void) function__;  // suppress unused var warning

    (void) function__;  // suppress unused var warning

    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    
    try {
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
      current_statement__ = 1;
      validate_non_negative_index("theta", "4", 4);
      std::vector<double> theta;
      theta = std::vector<double>(4, 0);
      
      current_statement__ = 2;
      validate_non_negative_index("z_init", "2", 2);
      std::vector<double> z_init;
      z_init = std::vector<double>(2, 0);
      
      current_statement__ = 3;
      validate_non_negative_index("sigma", "2", 2);
      std::vector<double> sigma;
      sigma = std::vector<double>(2, 0);
      
      current_statement__ = 4;
      validate_non_negative_index("z", "N", N);
      current_statement__ = 4;
      validate_non_negative_index("z", "2", 2);
      std::vector<std::vector<double>> z;
      z = std::vector<std::vector<double>>(N, std::vector<double>(2, 0));
      
      current_statement__ = 4;
      assign(z, nil_index_list(),
        integrate_ode_rk45(dz_dt_functor__(), z_init, 0, ts, theta,
          rep_array(0.0, 0), rep_array(0, 0), pstream__, 1e-5, 1e-3, 5e2),
        "assigning variable z");
      current_statement__ = 5;
      validate_non_negative_index("y_init", "2", 2);
      std::vector<double> y_init;
      y_init = std::vector<double>(2, 0);
      
      current_statement__ = 6;
      validate_non_negative_index("y", "N", N);
      current_statement__ = 6;
      validate_non_negative_index("y", "2", 2);
      std::vector<std::vector<double>> y;
      y = std::vector<std::vector<double>>(N, std::vector<double>(2, 0));
      
      current_statement__ = 7;
      assign(theta, cons_list(index_uni(1), nil_index_list()),
        normal_rng(1, 0.5, base_rng__), "assigning variable theta");
      current_statement__ = 8;
      assign(theta, cons_list(index_uni(3), nil_index_list()),
        normal_rng(1, 0.5, base_rng__), "assigning variable theta");
      current_statement__ = 9;
      assign(theta, cons_list(index_uni(2), nil_index_list()),
        normal_rng(0.05, 0.05, base_rng__), "assigning variable theta");
      current_statement__ = 10;
      assign(theta, cons_list(index_uni(4), nil_index_list()),
        normal_rng(0.05, 0.05, base_rng__), "assigning variable theta");
      current_statement__ = 11;
      assign(sigma, cons_list(index_uni(1), nil_index_list()),
        lognormal_rng(-1, 1, base_rng__), "assigning variable sigma");
      current_statement__ = 12;
      assign(sigma, cons_list(index_uni(2), nil_index_list()),
        lognormal_rng(-1, 1, base_rng__), "assigning variable sigma");
      current_statement__ = 13;
      assign(z_init, cons_list(index_uni(1), nil_index_list()),
        lognormal_rng(stan::math::log(10), 1, base_rng__),
        "assigning variable z_init");
      current_statement__ = 14;
      assign(z_init, cons_list(index_uni(2), nil_index_list()),
        lognormal_rng(stan::math::log(10), 1, base_rng__),
        "assigning variable z_init");
      current_statement__ = 18;
      for (size_t k = 1; k <= 2; ++k) {
        current_statement__ = 15;
        assign(y_init, cons_list(index_uni(k), nil_index_list()),
          lognormal_rng(stan::math::log(z_init[(k - 1)]), sigma[(k - 1)],
            base_rng__), "assigning variable y_init");
        current_statement__ = 16;
        assign(y,
          cons_list(index_omni(), cons_list(index_uni(k), nil_index_list())),
          lognormal_rng(
            stan::math::log(
              rvalue(z,
                cons_list(index_omni(),
                  cons_list(index_uni(k), nil_index_list())), "z")),
            sigma[(k - 1)], base_rng__), "assigning variable y");}
      current_statement__ = 1;
      for (size_t sym1__ = 1; sym1__ <= 4; ++sym1__) {
        current_statement__ = 1;
        current_statement__ = 1;
        check_greater_or_equal(function__, "theta[sym1__]",
                               theta[(sym1__ - 1)], 0);}
      current_statement__ = 2;
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        current_statement__ = 2;
        current_statement__ = 2;
        check_greater_or_equal(function__, "z_init[sym1__]",
                               z_init[(sym1__ - 1)], 0);}
      current_statement__ = 3;
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        current_statement__ = 3;
        current_statement__ = 3;
        check_greater_or_equal(function__, "sigma[sym1__]",
                               sigma[(sym1__ - 1)], 0);}
      current_statement__ = 6;
      for (size_t sym1__ = 1; sym1__ <= N; ++sym1__) {
        current_statement__ = 6;
        for (size_t sym2__ = 1; sym2__ <= 2; ++sym2__) {
          current_statement__ = 6;
          current_statement__ = 6;
          check_greater_or_equal(function__, "y[sym1__, sym2__]",
                                 y[(sym1__ - 1)][(sym2__ - 1)], 0);}}
      for (size_t sym1__ = 1; sym1__ <= 4; ++sym1__) {
        vars__.push_back(theta[(sym1__ - 1)]);}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        vars__.push_back(z_init[(sym1__ - 1)]);}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        vars__.push_back(sigma[(sym1__ - 1)]);}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        for (size_t sym2__ = 1; sym2__ <= N; ++sym2__) {
          vars__.push_back(z[(sym2__ - 1)][(sym1__ - 1)]);}}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        vars__.push_back(y_init[(sym1__ - 1)]);}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        for (size_t sym2__ = 1; sym2__ <= N; ++sym2__) {
          vars__.push_back(y[(sym2__ - 1)][(sym1__ - 1)]);}}
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array() 
    
  void transform_inits(const stan::io::var_context& context__,
                       std::vector<int>& params_i__,
                       std::vector<double>& vars__, std::ostream* pstream__) const {
    typedef double local_scalar_t__;
    vars__.resize(0);
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      
      pos__ = 1;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits() 
    
  void get_param_names(std::vector<std::string>& names__) const {
    
    names__.resize(0);
    names__.push_back("theta");
    names__.push_back("z_init");
    names__.push_back("sigma");
    names__.push_back("z");
    names__.push_back("y_init");
    names__.push_back("y");
    } // get_param_names() 
    
  void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.resize(0);
    std::vector<size_t> dims__;
    dims__.push_back(4);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dims__.push_back(2);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dims__.push_back(2);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dims__.push_back(N);
    
    dims__.push_back(2);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dims__.push_back(2);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dims__.push_back(N);
    
    dims__.push_back(2);
    dimss__.push_back(dims__);
    dims__.resize(0);
    
    } // get_dims() 
    
  void constrained_param_names(std::vector<std::string>& param_names__,
                               bool emit_transformed_parameters__ = true,
                               bool emit_generated_quantities__ = true) const {
    
    
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      for (size_t sym1__ = 1; sym1__ <= 4; ++sym1__) {
        {
          param_names__.push_back(std::string() + "theta" + '.' + std::to_string(sym1__));
        }}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.push_back(std::string() + "z_init" + '.' + std::to_string(sym1__));
        }}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.push_back(std::string() + "sigma" + '.' + std::to_string(sym1__));
        }}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          for (size_t sym2__ = 1; sym2__ <= N; ++sym2__) {
            {
              param_names__.push_back(std::string() + "z" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
            }}
        }}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.push_back(std::string() + "y_init" + '.' + std::to_string(sym1__));
        }}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          for (size_t sym2__ = 1; sym2__ <= N; ++sym2__) {
            {
              param_names__.push_back(std::string() + "y" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
            }}
        }}
    }
    
    } // constrained_param_names() 
    
  void unconstrained_param_names(std::vector<std::string>& param_names__,
                                 bool emit_transformed_parameters__ = true,
                                 bool emit_generated_quantities__ = true) const {
    
    
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      for (size_t sym1__ = 1; sym1__ <= 4; ++sym1__) {
        {
          param_names__.push_back(std::string() + "theta" + '.' + std::to_string(sym1__));
        }}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.push_back(std::string() + "z_init" + '.' + std::to_string(sym1__));
        }}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.push_back(std::string() + "sigma" + '.' + std::to_string(sym1__));
        }}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          for (size_t sym2__ = 1; sym2__ <= N; ++sym2__) {
            {
              param_names__.push_back(std::string() + "z" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
            }}
        }}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.push_back(std::string() + "y_init" + '.' + std::to_string(sym1__));
        }}
      for (size_t sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          for (size_t sym2__ = 1; sym2__ <= N; ++sym2__) {
            {
              param_names__.push_back(std::string() + "y" + '.' + std::to_string(sym2__) + '.' + std::to_string(sym1__));
            }}
        }}
    }
    
    } // unconstrained_param_names() 
    
  std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"theta\",\"type\":{\"name\":\"array\",\"length\":" << 4 << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"},{\"name\":\"z_init\",\"type\":{\"name\":\"array\",\"length\":" << 2 << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"},{\"name\":\"sigma\",\"type\":{\"name\":\"array\",\"length\":" << 2 << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"},{\"name\":\"z\",\"type\":{\"name\":\"array\",\"length\":" << N << ",\"element_type\":{\"name\":\"array\",\"length\":" << 2 << ",\"element_type\":{\"name\":\"real\"}}},\"block\":\"generated_quantities\"},{\"name\":\"y_init\",\"type\":{\"name\":\"array\",\"length\":" << 2 << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"},{\"name\":\"y\",\"type\":{\"name\":\"array\",\"length\":" << N << ",\"element_type\":{\"name\":\"array\",\"length\":" << 2 << ",\"element_type\":{\"name\":\"real\"}}},\"block\":\"generated_quantities\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"theta\",\"type\":{\"name\":\"array\",\"length\":" << 4 << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"},{\"name\":\"z_init\",\"type\":{\"name\":\"array\",\"length\":" << 2 << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"},{\"name\":\"sigma\",\"type\":{\"name\":\"array\",\"length\":" << 2 << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"},{\"name\":\"z\",\"type\":{\"name\":\"array\",\"length\":" << N << ",\"element_type\":{\"name\":\"array\",\"length\":" << 2 << ",\"element_type\":{\"name\":\"real\"}}},\"block\":\"generated_quantities\"},{\"name\":\"y_init\",\"type\":{\"name\":\"array\",\"length\":" << 2 << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"},{\"name\":\"y\",\"type\":{\"name\":\"array\",\"length\":" << N << ",\"element_type\":{\"name\":\"array\",\"length\":" << 2 << ",\"element_type\":{\"name\":\"real\"}}},\"block\":\"generated_quantities\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool emit_transformed_parameters__ = true,
                     bool emit_generated_quantities__ = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng__, params_r_vec, params_i_vec, vars_vec,
          emit_transformed_parameters__, emit_generated_quantities__, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    template <bool propto__, bool jacobian__, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto__,jacobian__,T_>(vec_params_r, vec_params_i, pstream);
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }

};
}

typedef dgp_lotka_volterra_model_namespace::dgp_lotka_volterra_model stan_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

#endif


