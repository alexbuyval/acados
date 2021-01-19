/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

#ifndef ACADOS_SOLVER_{{ model.name }}_H_
#define ACADOS_SOLVER_{{ model.name }}_H_

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

// ** global data **
// acados objects
struct all_global_data
{
    ocp_nlp_in * nlp_in;
    ocp_nlp_out * nlp_out;
    ocp_nlp_solver * nlp_solver;
    void * nlp_opts;
    ocp_nlp_plan * nlp_solver_plan;
    ocp_nlp_config * nlp_config;
    ocp_nlp_dims * nlp_dims;

// number of expected runtime parameters
    unsigned int nlp_np;

/* external functions */
// dynamics
    {% if solver_options.integrator_type == "ERK" %}
    external_function_param_casadi * forw_vde_casadi;
    external_function_param_casadi * expl_ode_fun;
    {% if solver_options.hessian_approx == "EXACT" %}
    external_function_param_casadi * hess_vde_casadi;
    {%- endif %}
    {% elif solver_options.integrator_type == "IRK" %}
    external_function_param_casadi * impl_dae_fun;
    external_function_param_casadi * impl_dae_fun_jac_x_xdot_z;
    external_function_param_casadi * impl_dae_jac_x_xdot_u_z;
    {%- if solver_options.hessian_approx == "EXACT" %}
    external_function_param_casadi * impl_dae_hess;
    {%- endif %}
    {% elif solver_options.integrator_type == "GNSF" %}
    external_function_param_casadi * gnsf_phi_fun;
    external_function_param_casadi * gnsf_phi_fun_jac_y;
    external_function_param_casadi * gnsf_phi_jac_y_uhat;
    external_function_param_casadi * gnsf_f_lo_jac_x1_x1dot_u_z;
    external_function_param_casadi * gnsf_get_matrices_fun;
    {% elif solver_options.integrator_type == "DISCRETE" %}
    external_function_param_casadi * discr_dyn_phi_fun;
    external_function_param_casadi * discr_dyn_phi_fun_jac_ut_xt;
    {%- if solver_options.hessian_approx == "EXACT" %}
    external_function_param_casadi * discr_dyn_phi_fun_jac_ut_xt_hess;
    {%- endif %}
    {%- endif %}

// cost
    {% if cost.cost_type == "NONLINEAR_LS" %}
    external_function_param_casadi * cost_y_fun;
    external_function_param_casadi * cost_y_fun_jac_ut_xt;
    external_function_param_casadi * cost_y_hess;
    {%- elif cost.cost_type == "EXTERNAL" %}
    external_function_param_casadi * ext_cost_fun;
    external_function_param_casadi * ext_cost_fun_jac;
    external_function_param_casadi * ext_cost_fun_jac_hess;
    {% endif %}
    {% if cost.cost_type_e == "NONLINEAR_LS" %}
    external_function_param_casadi cost_y_e_fun;
    external_function_param_casadi cost_y_e_fun_jac_ut_xt;
    external_function_param_casadi cost_y_e_hess;
    {% elif cost.cost_type_e == "EXTERNAL" %}
    external_function_param_casadi ext_cost_e_fun;
    external_function_param_casadi ext_cost_e_fun_jac;
    external_function_param_casadi ext_cost_e_fun_jac_hess;
    {%- endif %}

// constraints
    {%- if constraints.constr_type == "BGP" %}
    external_function_param_casadi * phi_constraint;
    {% elif constraints.constr_type == "BGH" and dims.nh > 0 %}
    external_function_param_casadi * nl_constr_h_fun_jac;
    external_function_param_casadi * nl_constr_h_fun;
    external_function_param_casadi * nl_constr_h_fun_jac_hess;
    {% endif %}

    {% if constraints.constr_type_e == "BGP" %}
    external_function_param_casadi phi_e_constraint;
    {% elif constraints.constr_type_e == "BGH" and dims.nh_e > 0 %}
    external_function_param_casadi nl_constr_h_e_fun_jac;
    external_function_param_casadi nl_constr_h_e_fun;
    external_function_param_casadi nl_constr_h_e_fun_jac_hess;
    {%- endif %}
};

#ifdef __cplusplus
extern "C" {
#endif

int acados_create(int N, double time_step, struct all_global_data* data);
int acados_update_params(int stage, double *value, int np, struct all_global_data* data);
int acados_solve(struct all_global_data* data);
int acados_free(int N, struct all_global_data* data);
void acados_print_stats(struct all_global_data* data);

ocp_nlp_in * acados_get_nlp_in(struct all_global_data* data);
ocp_nlp_out * acados_get_nlp_out(struct all_global_data* data);
ocp_nlp_solver * acados_get_nlp_solver(struct all_global_data* data);
ocp_nlp_config * acados_get_nlp_config(struct all_global_data* data);
void * acados_get_nlp_opts(struct all_global_data* data);
ocp_nlp_dims * acados_get_nlp_dims(struct all_global_data* data);
ocp_nlp_plan * acados_get_nlp_plan(struct all_global_data* data);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SOLVER_{{ model.name }}_H_
