use osqp_rust_sys as ffi;
use std::fmt;
use std::slice;
use std::time::Duration;

use crate::{float, Problem};

/// The result of solving a problem.
#[derive(Clone, Debug)]
pub enum Status<'a> {
    Solved(Solution<'a>),
    SolvedInaccurate(Solution<'a>),
    MaxIterationsReached(Solution<'a>),
    TimeLimitReached(Solution<'a>),
    PrimalInfeasible(PrimalInfeasibilityCertificate<'a>),
    PrimalInfeasibleInaccurate(PrimalInfeasibilityCertificate<'a>),
    DualInfeasible(DualInfeasibilityCertificate<'a>),
    DualInfeasibleInaccurate(DualInfeasibilityCertificate<'a>),
    NonConvex(Failure<'a>),
    // Prevent exhaustive enum matching
    #[doc(hidden)]
    __Nonexhaustive,
}

/// A solution to a problem.
#[derive(Clone)]
pub struct Solution<'a> {
    prob: &'a Problem,
}

/// A proof of primal infeasibility.
#[derive(Clone)]
pub struct PrimalInfeasibilityCertificate<'a> {
    prob: &'a Problem,
}

/// A proof of dual infeasibility.
#[derive(Clone)]
pub struct DualInfeasibilityCertificate<'a> {
    prob: &'a Problem,
}

/// A problem that failed to be solved.
#[derive(Clone)]
pub struct Failure<'a> {
    prob: &'a Problem,
}

/// The status of the polish operation.
#[derive(Copy, Clone, Debug, Hash, PartialEq)]
pub enum PolishStatus {
    Successful,
    Unsuccessful,
    Unperformed,
    // Prevent exhaustive enum matching
    #[doc(hidden)]
    __Nonexhaustive,
}

impl<'a> Status<'a> {
    pub(crate) fn from_problem(prob: &'a Problem) -> Status<'a> {
        use std::os::raw::c_int;
        unsafe {
            match (*(*prob.workspace).info).status_val as c_int {
                ffi::src::src::auxil::OSQP_SOLVED => Status::Solved(Solution { prob }),
                ffi::src::src::auxil::OSQP_SOLVED_INACCURATE => Status::SolvedInaccurate(Solution { prob }),
                ffi::src::src::auxil::OSQP_MAX_ITER_REACHED => Status::MaxIterationsReached(Solution { prob }),
                ffi::src::src::auxil::OSQP_TIME_LIMIT_REACHED => Status::TimeLimitReached(Solution { prob }),
                ffi::src::src::auxil::OSQP_PRIMAL_INFEASIBLE => {
                    Status::PrimalInfeasible(PrimalInfeasibilityCertificate { prob })
                }
                ffi::src::src::auxil::OSQP_PRIMAL_INFEASIBLE_INACCURATE => {
                    Status::PrimalInfeasibleInaccurate(PrimalInfeasibilityCertificate { prob })
                }
                ffi::src::src::auxil::OSQP_DUAL_INFEASIBLE => {
                    Status::DualInfeasible(DualInfeasibilityCertificate { prob })
                }
                ffi::src::src::auxil::OSQP_DUAL_INFEASIBLE_INACCURATE => {
                    Status::DualInfeasibleInaccurate(DualInfeasibilityCertificate { prob })
                }
                ffi::src::src::auxil::OSQP_NON_CVX => Status::NonConvex(Failure { prob }),
                _ => unreachable!(),
            }
        }
    }

    /// Returns the primal variables at the solution if the problem is `Solved`.
    pub fn x(&self) -> Option<&'a [float]> {
        self.solution().map(|s| s.x())
    }

    /// Returns the solution if the problem is `Solved`.
    pub fn solution(&self) -> Option<Solution<'a>> {
        match *self {
            Status::Solved(ref solution) => Some(solution.clone()),
            _ => None,
        }
    }

    /// Returns the number of iterations taken by the solver.
    pub fn iter(&self) -> u32 {
        unsafe {
            // cast safe as more than 2 billion iterations would be unreasonable
            (*(*self.prob().workspace).info).iter as u32
        }
    }

    /// Returns the time taken for the setup phase.
    pub fn setup_time(&self) -> Duration {
        unsafe { secs_to_duration((*(*self.prob().workspace).info).setup_time) }
    }

    /// Returns the time taken for the solve phase.
    pub fn solve_time(&self) -> Duration {
        unsafe { secs_to_duration((*(*self.prob().workspace).info).solve_time) }
    }

    /// Returns the time taken for the polish phase.
    pub fn polish_time(&self) -> Duration {
        unsafe { secs_to_duration((*(*self.prob().workspace).info).polish_time) }
    }

    /// Returns the total time taken by the solver.
    ///
    /// This includes the time taken for the setup phase on the first solve.
    pub fn run_time(&self) -> Duration {
        unsafe { secs_to_duration((*(*self.prob().workspace).info).run_time) }
    }

    /// Returns the number of rho updates.
    pub fn rho_updates(&self) -> u32 {
        unsafe {
            // cast safe as more than 2 billion updates would be unreasonable
            (*(*self.prob().workspace).info).rho_updates as u32
        }
    }

    /// Returns the current best estimate of rho.
    pub fn rho_estimate(&self) -> float {
        unsafe { (*(*self.prob().workspace).info).rho_estimate }
    }

    fn prob(&self) -> &'a Problem {
        match *self {
            Status::Solved(ref solution)
            | Status::SolvedInaccurate(ref solution)
            | Status::MaxIterationsReached(ref solution)
            | Status::TimeLimitReached(ref solution) => solution.prob,
            Status::PrimalInfeasible(ref cert) | Status::PrimalInfeasibleInaccurate(ref cert) => {
                cert.prob
            }
            Status::DualInfeasible(ref cert) | Status::DualInfeasibleInaccurate(ref cert) => {
                cert.prob
            }
            Status::NonConvex(ref failure) => failure.prob,
            Status::__Nonexhaustive => unreachable!(),
        }
    }
}

impl<'a> Solution<'a> {
    /// Returns the primal variables at the solution.
    pub fn x(&self) -> &'a [float] {
        unsafe { slice::from_raw_parts((*(*self.prob.workspace).solution).x, self.prob.n) }
    }

    /// Returns the dual variables at the solution.
    ///
    /// These are the Lagrange multipliers of the constraints `l <= Ax <= u`.
    pub fn y(&self) -> &'a [float] {
        unsafe { slice::from_raw_parts((*(*self.prob.workspace).solution).y, self.prob.m) }
    }

    /// Returns the status of the polish operation.
    pub fn polish_status(&self) -> PolishStatus {
        unsafe {
            match (*(*self.prob.workspace).info).status_polish {
                1 => PolishStatus::Successful,
                -1 => PolishStatus::Unsuccessful,
                0 => PolishStatus::Unperformed,
                _ => unreachable!(),
            }
        }
    }

    /// Returns the primal objective value.
    pub fn obj_val(&self) -> float {
        unsafe { (*(*self.prob.workspace).info).obj_val }
    }

    /// Returns the norm of primal residual.
    pub fn pri_res(&self) -> float {
        unsafe { (*(*self.prob.workspace).info).pri_res }
    }

    /// Returns the norm of dual residual.
    pub fn dua_res(&self) -> float {
        unsafe { (*(*self.prob.workspace).info).dua_res }
    }
}

impl<'a> fmt::Debug for Solution<'a> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("Solution")
            .field("x", &self.x())
            .field("y", &self.y())
            .field("polish_status", &self.polish_status())
            .field("obj_val", &self.obj_val())
            .field("pri_res", &self.pri_res())
            .field("dua_res", &self.dua_res())
            .finish()
    }
}

impl<'a> PrimalInfeasibilityCertificate<'a> {
    /// Returns the certificate of primal infeasibility.
    ///
    /// For further explanation see [Infeasibility detection in the alternating direction method of
    /// multipliers for convex
    /// optimization](http://www.optimization-online.org/DB_HTML/2017/06/6058.html).
    pub fn delta_y(&self) -> &'a [float] {
        unsafe { slice::from_raw_parts((*self.prob.workspace).delta_y, self.prob.m) }
    }
}

impl<'a> fmt::Debug for PrimalInfeasibilityCertificate<'a> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("PrimalInfeasibilityCertificate")
            .field("delta_y", &self.delta_y())
            .finish()
    }
}

impl<'a> DualInfeasibilityCertificate<'a> {
    /// Returns the certificate of dual infeasibility.
    ///
    /// For further explanation see [Infeasibility detection in the alternating direction method of
    /// multipliers for convex
    /// optimization](http://www.optimization-online.org/DB_HTML/2017/06/6058.html).
    pub fn delta_x(&self) -> &'a [float] {
        unsafe { slice::from_raw_parts((*self.prob.workspace).delta_x, self.prob.n) }
    }
}

impl<'a> fmt::Debug for DualInfeasibilityCertificate<'a> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("DualInfeasibilityCertificate")
            .field("delta_x", &self.delta_x())
            .finish()
    }
}

impl<'a> fmt::Debug for Failure<'a> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("Failure").finish()
    }
}

fn secs_to_duration(secs: float) -> Duration {
    let whole_secs = secs.floor() as u64;
    let nanos = (secs.fract() * 1e9) as u32;
    Duration::new(whole_secs, nanos)
}
