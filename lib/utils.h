#ifndef UTILS
#define UTILS

#include "vectors.h"

#include <fstream>

// TODO: Add documentation
void log_units();

void log_progress(int current_step, int tot_steps);

void save_pos_state(std::ofstream *file, double t, Vec3* xs, int n_part);


#endif
