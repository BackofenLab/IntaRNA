#ifndef _ENERGY__
#define _ENERGY__

#include <stdlib.h>
#include <iostream>
#include "basics.h"

using namespace std;

void   scale_parameters();
void   scaling_function(double scaling_factor, double& energy, double enthalpy);
double calc_ensemble_free_energy(char* sequence, int start_unfold, int end_unfold);
double calc_ED (char* sequence, int start_unpaired, int end_unpaired);
double calc_PU (double ED);
void   calc_ED_target_matrix();
void   calc_ED_miRNA_matrix();
void calc_ED_miRNA_matrix_up();
void   calc_ED_target_matrix_window();
void   calc_ED_target_matrix_plfold(char* sequence);
//void   calc_ED_target_matrix_old_plfold(char* sequence);
double AU_penalty(char,char);
// double ed5(char, char, char); 
// double ed3(char, char, char);

//bool StackEnd(int bp_pos);
//double Zero_or_StemEndAU(int bp_pos, int bp_assign);

#endif   // _ENERGY_
