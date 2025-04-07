#ifndef EXTRACTOR_H
#define EXTRACTOR_H

// Variable definitions
#include "../Vars/simparam.h"

bool extract_boolean_value(const char *json, const char *key, bool *value);
char *read_file(const char *filename);
int extract_magnetic_field(const char *json);
int extract_toroidal_machine_geometry(const char *json);
int extract_neutral_pressure(const char *json);
int extract_rf_power(const char *json);
int extract_type(const char *json);
int extract_general_ec(const char *json);
int extract_necfix(const char *json);
int extract_tomas(const char *json);
int extract_general_ic(const char *json);
int extract_at_lhr(const char *json);
int extract_other(const char *json);
int extract_diffusion(const char *json);
int extract_convection(const char *json);
int extract_tune_transport(const char *json);
int extract_physics_to_include(const char *json);
int extract_initial_conditions(const char *json);
int extract_edge_conditions(const char *json);
int extract_simulation_grid(const char *json);
int extract_input_file(const char *json);
int extract_time_step(const char *json);
int extract_time_step_for_rf_coupling(const char *json);
int extract_advanced_time_step_settings(const char *json);
int extract_output_parameters(const char *json);
int extract_solver_parameters(const char *json);
int extract_input_file_name(const char *json);
int extract_bmanuel_input(const char *json);

#endif // EXTRACTOR_H
