use osqp_rust_sys as ffi;
use std::mem;
use std::ptr;
use std::time::Duration;

use crate::{float, Problem};

/// The linear system solver for OSQP to use.
#[derive(Clone, Debug, PartialEq)]
pub enum LinsysSolver {
    Qdldl,
    MklPardiso,
    // Prevent exhaustive enum matching
    #[doc(hidden)]
    __Nonexhaustive,
}

macro_rules! u32_to_osqp_int {
    ($name:ident, $value:expr) => {{
        let value = $value;
        assert!(
            value as u64 <= ffi::src::src::osqp::c_int::max_value() as u64,
            "{} must be smaller than the largest isize value",
            stringify!($name)
        );
        value as ffi::src::src::osqp::c_int
    }};
}

macro_rules! rust_type {
    (float) => (float);
    (u32) => (u32);
    (option_u32) => (Option<u32>);
    (bool) => (bool);
    (linsys_solver) => (LinsysSolver);
    (option_duration) => (Option<Duration>);
}

macro_rules! convert_rust_type {
    ($name:ident, float, $value:expr) => ($value);
    ($name:ident, u32, $value:expr) => (u32_to_osqp_int!($name, $value));
    ($name:ident, option_u32, $value:expr) => (u32_to_osqp_int!($name, $value.unwrap_or(0)));
    ($name:ident, bool, $value:expr) => ($value as ffi::src::src::osqp::c_int);
    ($name:ident, linsys_solver, $value:expr) => (
        match $value {
            LinsysSolver::Qdldl => ffi::src::src::osqp::QDLDL_SOLVER,
            LinsysSolver::MklPardiso => ffi::src::src::osqp::MKL_PARDISO_SOLVER,
            LinsysSolver::__Nonexhaustive => unreachable!(),
        }
    );
    ($name:ident, option_duration, $value:expr) => (
        $value.map(|v| {
            let mut secs = duration_to_secs(v);
            // Setting time_limit to 0.0 disables the time limit to we treat a duration of zero as
            // a very small time limit instead.
            if secs == 0.0 {
                secs = 1e-12;
            }
            secs
        }).unwrap_or(0.0)
    );
}

macro_rules! settings {
    ($problem_ty:ty, $(
        #[$doc:meta] $name:ident: $typ:ident $([$update_name:ident, $update_ffi:ident])*,
    )*) => (
        /// The settings used when initialising a solver.
        pub struct Settings {
            pub(crate) inner: ffi::src::src::osqp::OSQPSettings,
        }

        impl Settings {
            $(
                #[$doc]
                pub fn $name(mut self, value: rust_type!($typ)) -> Settings {
                    self.inner.$name = convert_rust_type!($name, $typ, value);
                    Settings {
                        inner: self.inner
                    }
                }
            )*
        }

        impl Clone for Settings {
            fn clone(&self) -> Settings {
                unsafe {
                    Settings {
                        inner: ptr::read(&self.inner)
                    }
                }
            }
        }

        impl Default for Settings {
            fn default() -> Settings {
                unsafe {
                    let mut settings: ffi::src::src::osqp::OSQPSettings = mem::zeroed();
                    ffi::src::src::osqp::osqp_set_default_settings(&mut settings);
                    Settings {
                        inner: settings
                    }
                }
            }
        }

        unsafe impl Send for Settings {}
        unsafe impl Sync for Settings {}

        impl $problem_ty {
            $($(
                #[$doc]
                pub fn $update_name(&mut self, value: rust_type!($typ)) {
                    unsafe {
                        let ret = ffi::src::src::osqp::$update_ffi(
                            self.workspace,
                            convert_rust_type!($name, $typ, value)
                        );
                        if ret != 0 {
                            panic!("updating {} failed", stringify!($name));
                        }
                    }
                }
            )*)*
        }
    );
}

settings! {
    Problem,

    #[doc = "Sets the ADMM step rho."]
    rho: float [update_rho, osqp_update_rho],

    #[doc = "Sets the ADMM step sigma."]
    sigma: float,

    #[doc = "
    Sets the number of heuristic data scaling iterations.

    If `None` scaling is disabled.

    Panics on 32-bit platforms if the value is above `i32::max_value()`.
    "]
    scaling: option_u32,

    #[doc = "Enables choosing rho adaptively."]
    adaptive_rho: bool,

    #[doc = "
    Sets the number of iterations between rho adaptations.

    If `None` it is automatic.

    Panics on 32-bit platforms if the value is above `i32::max_value()`.
    "]
    adaptive_rho_interval: option_u32,

    #[doc = "
    Sets the tolerance for adapting rho.

    The new rho has to be `value` times larger or `1/value` times smaller than the current rho to
    trigger a new factorization.
    "]
    adaptive_rho_tolerance: float,

    #[doc = "Set the interval for adapting rho as a fraction of the setup time."]
    adaptive_rho_fraction: float,

    #[doc = "
    Sets the maximum number of ADMM iterations.

    Panics on 32-bit platforms if the value is above `i32::max_value()`.
    "]
    max_iter: u32 [update_max_iter, osqp_update_max_iter],

    #[doc = "Sets the absolute convergence tolerance."]
    eps_abs: float [update_eps_abs, osqp_update_eps_abs],

    #[doc = "Sets the relative convergence tolerance."]
    eps_rel: float [update_eps_rel, osqp_update_eps_rel],

    #[doc = "Sets the primal infeasibility tolerance."]
    eps_prim_inf: float [update_eps_prim_inf, osqp_update_eps_prim_inf],

    #[doc = "Sets the dual infeasibility tolerance."]
    eps_dual_inf: float [update_eps_dual_inf, osqp_update_eps_dual_inf],

    #[doc = "Sets the linear solver relaxation parameter."]
    alpha: float [update_alpha, osqp_update_alpha],

    #[doc = "Sets the linear system solver to use."]
    linsys_solver: linsys_solver,

    #[doc = "Sets the polishing regularization parameter."]
    delta: float [update_delta, osqp_update_delta],

    #[doc = "Enables polishing the ADMM solution."]
    polish: bool [update_polish, osqp_update_polish],

    #[doc = "
    Sets the number of iterative refinement steps to use when polishing.

    Panics on 32-bit platforms if the value is above `i32::max_value()`.
    "]
    polish_refine_iter: u32 [update_polish_refine_iter, osqp_update_polish_refine_iter],

    #[doc = "Enables writing progress to stdout."]
    verbose: bool [update_verbose, osqp_update_verbose],

    #[doc = "Enables scaled termination criteria."]
    scaled_termination: bool [update_scaled_termination, osqp_update_scaled_termination],

    #[doc = "
    Sets the number of ADMM iterations between termination checks.

    If `None` termination checking is disabled.

    Panics on 32-bit platforms if the value is above `i32::max_value()`.
    "]
    check_termination: option_u32 [update_check_termination, osqp_update_check_termination],

    #[doc = "Enables warm starting the primal and dual variables from the previous solution."]
    warm_start: bool [update_warm_start, osqp_update_warm_start],

    #[doc = "Sets the solve time limit."]
    time_limit: option_duration [update_time_limit, osqp_update_time_limit],
}

fn duration_to_secs(dur: Duration) -> float {
    dur.as_secs() as float + dur.subsec_nanos() as float * 1e-9
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(target_pointer_width = "32")]
    #[test]
    #[should_panic]
    fn large_u32_settings_value_panics_on_32_bit() {
        // Larger than i32::max_value()
        Settings::default().polish_refine_iter(3_000_000_000);
    }

    #[test]
    fn duration_to_secs_examples() {
        assert_eq!(duration_to_secs(Duration::new(2, 0)), 2.0);
        assert_eq!(duration_to_secs(Duration::new(8, 100_000_000)), 8.1);
        assert_eq!(duration_to_secs(Duration::new(0, 10_000_000)), 0.01);
    }
}
