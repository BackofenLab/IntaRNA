#include "energy.h"

#include "./data/stacking_energies.dat"
#include "./data/stacking_enthalpies.dat"
#include "./data/interior_loop_1_1_energies.dat"
#include "./data/interior_loop_1_1_enthalpies.dat"
#include "./data/interior_loop_1_2_energies.dat"
#include "./data/interior_loop_1_2_enthalpies.dat"
#include "./data/interior_loop_2_2_energies.dat"
#include "./data/interior_loop_2_2_enthalpies.dat"
#include "./data/loop_destabilizing_energies.dat"
#include "./data/single_base_stacking_energies.dat"
#include "./data/single_base_stacking_enthalpies.dat"
#include "./data/terminal_mismatch_hairpin.dat"
#include "./data/terminal_mismatch_interior.dat"
#include "./data/terminal_mismatch_enthalpies.dat"

void pup_factor_correct(double *pup,int length,int unpaired_length); //adapted of a Vienna function


/**********************************************************************************
*        scales all free energy parameter for a temperature different from 37°C   *
**********************************************************************************/

void scale_parameters()
{
    double scaling_factor = (temperature+K_0)/T_MEASURE; // temperature scaling factor

    // scale free energy parameter
    for(size_t i = 0; i<sizeof(stacking_energies)/sizeof(double); i++) {
        scaling_function(scaling_factor,stacking_energies[i],stacking_enthalpies[i]);
    }
    for(size_t i = 0; i<sizeof(interior_loop_1_1_energy)/sizeof(double); i++) {
        scaling_function(scaling_factor,interior_loop_1_1_energy[i],interior_loop_1_1_enthalpy[i]);
    }
    for(size_t i = 0; i<sizeof(interior_loop_1_2_energy)/sizeof(double); i++) {
        scaling_function(scaling_factor,interior_loop_1_2_energy[i],interior_loop_1_2_enthalpy[i]);
    }
    for(size_t i = 0; i<sizeof(interior_loop_2_2_energy)/sizeof(double); i++) {
        scaling_function(scaling_factor,interior_loop_2_2_energy[i],interior_loop_2_2_enthalpy[i]);
    }
    for(size_t i = 0; i<sizeof(single_base_stacking_energy)/sizeof(double); i++) {
        scaling_function(scaling_factor,single_base_stacking_energy[i],single_base_stacking_enthalpies[i]);
        // dangling end contributions must be <= 0
        if (single_base_stacking_energy[i] > 0) {
            single_base_stacking_energy[i] = 0;
        }
    }
    for(size_t i = 0; i<sizeof(mismatch_energies_hairpin)/sizeof(double); i++) {
        scaling_function(scaling_factor,mismatch_energies_hairpin[i],mismatch_enthalpies[i]);
    }
    for(size_t i = 0; i<sizeof(mismatch_energies_interior)/sizeof(double); i++) {
        scaling_function(scaling_factor,mismatch_energies_interior[i],mismatch_enthalpies[i]);
    }
    for(size_t i = 0; i<sizeof(loop_destabilizing_energies)/sizeof(double); i++) {
        loop_destabilizing_energies[i] *= scaling_factor;
    }

    loop_extrapolation *= scaling_factor;
    ninio_correction *= scaling_factor;
    duplex_init *= scaling_factor;
}


/**********************************************************************************
*        scales passed energy parameter to given temperature (in [K])             *
**********************************************************************************/

void scaling_function(double scaling_factor, double& energy, double enthalpy)
{
    energy = enthalpy - (scaling_factor * (enthalpy - energy));
}


/**********************************************************************************
*         calulates ensemble energies (constraints possible)                      *
*        if unconstraint: start_unfold = -1 and end_unfold = -1                   *
**********************************************************************************/

double calc_ensemble_free_energy(char* sequence, int start_unfold, int end_unfold)
{
   int len = (int)strlen(sequence);
//    double kT = (temperature+273.15)*1.98717/1000.0; /* in Kcal */
   double sfact = 1.07;

   //structure, might be constraint (if not: only '.'s, which means 'unconstraint')
   char structure[len+1];
   char c_structure[len+1];
   for (int i=0; i<len; i++)
   {
      //  we use positions starting at 0
      if ((start_unfold <= i) && (i <= end_unfold))
         c_structure[i] = 'x';
      else
         c_structure[i] = '.';
   }
   c_structure[len] = '\0';

   strncpy(structure, c_structure, len+1);

   double min_en = fold(sequence, structure);
   //   printf("Sequence:      %s\nConstraints:   %s\nmfe structure: %s\n", sequence, c_structure, structure);
   //   printf("\n minimum free energy = %6.2f kcal/mol\n\n", min_en);

   char pf_structure[len+1];
   if (dangles==1) 
   {
      dangles=2;   /* recompute with dangles as in pf_fold() */
      min_en = energy_of_struct(sequence, structure);
      dangles=1;
   }

   free_arrays ();

   pf_scale = exp(-(sfact*min_en)/RT/len);
   if (len>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);

   init_pf_fold(len);
   strncpy(pf_structure, c_structure, len+1);
   double energy = pf_fold(sequence, pf_structure);
   //   printf(" free energy of ensemble = %6.2f kcal/mol\n", energy);
   //   printf(" frequency of mfe structure in ensemble %g; \n\n", exp((energy-min_en)/kT));
   free_pf_arrays();

   return(energy);
}


/******************************************************************************
* calculates ED = E_unpaired(start_unpaired,end_unfold) - E_all)              *
******************************************************************************/

double calc_ED (char* sequence, int start_unpaired, int end_unpaired)
{
   double E_all = calc_ensemble_free_energy(sequence, -1, -1);
   double E_unpaired = calc_ensemble_free_energy(sequence, start_unpaired, end_unpaired);
   return E_unpaired-E_all;
} 

/******************************************************************************
* calculates ED values for all intervals in targetRNA                         *
******************************************************************************/

void calc_ED_target_matrix()
{
  double E_all = calc_ensemble_free_energy(sequence_target, -1, -1);

  //cout << "no window, no RNAplfold" << endl;
  string prefix =" ED(targetRNA) computation: " ;    
  int count=0;

  for(int i=0;i<(int)strlen(sequence_target);i++)
    {
    //for(int j=i;j<(int)strlen(sequence_target);j++)
    //Mininum is needed, if no -w option is given but a -l option
    for(int j=i; j < Minimum((int)strlen(sequence_target),i+max_unpaired_length); j++)
      {
	ED_target[i][j]= calc_ensemble_free_energy(sequence_target,i,j)-E_all;

        // Ouput progress only if no sequence file given
        if( target_file==false && miRNA_file==false)
         { 
           count++;
           //number of calculated values out of all values (times 200 since: *100 to change to percent and *2 since actually the denominator has to be divided by 2)
           //the denominator represents the number of values in a large upper triangle matrix minus the number in a smaller upper triangle matrix
           cout << prefix << count*200/(((int)strlen(sequence_target) * ((int)strlen(sequence_target) + 1))- (((int)strlen(sequence_target)-max_unpaired_length)*((int)strlen(sequence_target)-max_unpaired_length+1))) << "%\r" ;
	   cout.flush();
         }
      }
    }
  
  if( target_file==false && miRNA_file==false)
  {
     cout << endl;
  }
}

/******************************************************************************
* calculates ED values for all intervals in miRNA                             *
******************************************************************************/

void calc_ED_miRNA_matrix()
{
  //change direction of the second sequence again for a correct calculation of
  //the ED-values
  change_direction(sequence_miRNA);
  double E_all = calc_ensemble_free_energy(sequence_miRNA , -1, -1);

  int miRNA_len = (int)strlen(sequence_miRNA);
  string prefix =" ED(    ncRNA) computation: " ;    
  int count=0;
  for(int i=0; i<miRNA_len; i++)
    {
    for(int j=i; j<miRNA_len; j++)
      {
        //since we need the ED-values of the original sequence, but we use the inverse
        // one, we have to convert this
	ED_miRNA [miRNA_len-1-j][miRNA_len-1-i]= calc_ensemble_free_energy(sequence_miRNA,i,j)-E_all;
	//ED_miRNA [i][j]= calc_ensemble_free_energy(sequence_miRNA,i,j)-E_all;

        // Ouput progress only if no sequence file given
        if( target_file==false && miRNA_file==false)
         { 
           count++;
           cout << prefix << count*200/(miRNA_len*(miRNA_len+1)) << "%\r" ;
	   cout.flush();
         }
      }
    }

  if( target_file==false && miRNA_file==false)
    {
      cout << endl << endl;
    }
  //change direction back
  change_direction(sequence_miRNA);
} 

/******************************************************************************
* calculates ED values for all intervals <= max_unpaired_length of miRNA  *
* using RNAup                                                             *
******************************************************************************/

void calc_ED_miRNA_matrix_up()
{

	// structure with four arrays containing the contributions for PU of maximum length w
	// arrangement: i = 1..seq_length, j = 0..w-1 (length of unpaired region - 1)
	pu_contrib *pu_struct;

	double pu;
	double sfact = 1.07;

  //change direction of the second sequence again for a correct calculation of
  //the ED-values
  change_direction(sequence_miRNA);
  double E_all = calc_ensemble_free_energy(sequence_miRNA , -1, -1);

  int miRNA_len = (int)strlen(sequence_miRNA);

	// a scaling factor used by pf_fold() to avoid overflows
	pf_scale = exp(-(sfact*E_all)/RT/miRNA_len);
	// allocates memory for folding sequences not longer than length
	// sets up pairing matrix and energy parameters
	// has to be called before the first call to pf_fold()
	init_pf_fold(miRNA_len);

	// calculate partition function
	pf_fold(sequence_miRNA, NULL);
	// calculate unpaired probabilities for all intervals <= Minimum(miRNA_len,max_unpaired_length)
	pu_struct = pf_unstru(sequence_miRNA, Minimum(miRNA_len,max_unpaired_length));

  for(int u=1;u<=Minimum(miRNA_len,max_unpaired_length);u++) {
    for(int i=0; i<miRNA_len-u+1; i++) {
			// total probability of being unpaired is sum of 4 contributions
			// first index: position of unpaired region start counting from 1
			// second index: length of unpaired region - 1
			pu = pu_struct->H[i+1][u-1]+pu_struct->I[i+1][u-1]+pu_struct->M[i+1][u-1]+pu_struct->E[i+1][u-1];

			// since we need the ED-values of the original sequence, but we use the inverse
			// one, we have to convert this
			if (pu == 0.0)  // is not exact since we have double values, but pu can get really small
				ED_miRNA[miRNA_len-1-(i+u-1)][miRNA_len-1-i] = 
					(fabs(ED_miRNA_scaling - 0.0) < DOUBLE_DIFF) ? 0 : MAX_DOUBLE;
			else
				ED_miRNA[miRNA_len-1-(i+u-1)][miRNA_len-1-i] = -RT * ED_miRNA_scaling * log(pu); 

      }
    }

//   for(int i=0; i<miRNA_len; i++)
//     for(int j=i; j<miRNA_len; j++)
// 			cout << "ED_miRNA[" <<  i << "][" << j << "]: " << ED_miRNA[i][j] << endl;

	free_pu_contrib(pu_struct);
  // free space allocated for pf_fold()
	free_pf_arrays();

  //change direction back
  change_direction(sequence_miRNA);
} 

/******************************************************************************
* calculates ED values for all intervals in window of targetRNA               *
******************************************************************************/

void calc_ED_target_matrix_window()
{
  string prefix       = " ED(targetRNA) computation: " ; 
  int    count     = 0;

  //cout << "window without RNAplfold used" << endl;

  // compute left ED values
  char* subseq=(char*) malloc(sizeof(char)*(window_size+1));
  for(int k=0;k<window_size;k++)
    subseq[k]=sequence_target[k];
  subseq[window_size] = '\0';

  double E_all = calc_ensemble_free_energy(subseq, -1, -1);

  for(int d=1;d<=max_unpaired_length;d++)
    {
      for (int i=0;i<(window_size-d)/2;i++)
	{
	  int j=i+d-1;
	  ED_target[i][j]= calc_ensemble_free_energy(subseq,i,j)-E_all;

          // Ouput progress only if no sequence file given
          if( target_file==false && miRNA_file==false)
           { 
             count++;
             cout << prefix << count*66/(max_unpaired_length*window_size) << "%\r" ;
	     cout.flush();
           }
	}
    }
  free(subseq);

  // compute right ED values
  subseq=(char*) malloc(sizeof(char)*(window_size+1));
  for(int k=0;k<window_size;k++)
    subseq[k]=sequence_target[k+(strlen(sequence_target)-window_size)];
  subseq[window_size] = '\0';

  E_all = calc_ensemble_free_energy(subseq, -1, -1);
  
  for(int d=1;d<=max_unpaired_length;d++)
    {
      for (int i=(strlen(sequence_target)-window_size+1+(window_size-d)/2);i<(int)(strlen(sequence_target)-d+1);i++)
	{
	  int j=i+d-1;
	  
	  int i_in_window = i - (int)(strlen(sequence_target)-window_size+1+(window_size-d)/2);
	  int j_in_window = j - (int)(strlen(sequence_target)-window_size+1+(window_size-d)/2);
	  
	  ED_target[i][j]= calc_ensemble_free_energy(subseq,i_in_window,j_in_window)-E_all;

          // Ouput progress only if no sequence file given
          if( target_file==false && miRNA_file==false)
           { 
             count++;
             cout << prefix << count*66/(max_unpaired_length*window_size) << "%\r" ;
	     cout.flush();
           }
	}
    }
  free(subseq);

  // compute regular ED values
  int per_cent = count*66/(max_unpaired_length*window_size);
  count = 0;

  for (int w1=0;w1<(int)(strlen(sequence_target)-window_size+1);w1++)
    {
      int w2=w1+window_size-1;

      // RNA_target subsequence from w1 to w2
      subseq=(char*) malloc(sizeof(char)*(window_size+1));
      for(int k=w1;k<=w2;k++)
	subseq[k-w1]=sequence_target[k];
      subseq[window_size] = '\0';

      E_all = calc_ensemble_free_energy(subseq, -1, -1);

      for(int d=1;d<=max_unpaired_length;d++)
	{
	  int i=w1+(window_size-d)/2;
	  int j=i+d-1;
	  
	  int i_in_window = i - w1;
	  int j_in_window = j - w1;

	  ED_target[i][j]= calc_ensemble_free_energy(subseq,i_in_window,j_in_window)-E_all;
	  /*  cout << "subseq: " << subseq << endl;
	     cout << "Eunpaired_i_j      : " << calc_ensemble_free_energy(subseq, i, j) << endl;
	     cout << "Eunpaired_i_in_j_in: " << calc_ensemble_free_energy(subseq, i_in_window, j_in_window) << endl;
	     cout << "Eall               : " << E_all << endl;
	     cout << "i: " << i << " j: " << j << " ED_target[i][j]: " << ED_target[i][j] << endl;
	     cout << "i_in: " << i_in_window << " j_in: " << j_in_window << endl << endl;*/

          // Ouput progress only if no sequence file given
          if( target_file==false && miRNA_file==false)
           { 
             count++;
             cout << prefix << (count*(100-per_cent)/((int)(strlen(sequence_target)-window_size+1)*max_unpaired_length))+per_cent << "%\r" ;
	     cout.flush();
           }
	}
      free(subseq);
    }
  if( target_file==false && miRNA_file==false)
    {
      cout << endl;
    }
} 


/******************************************************************************
* calculates ED values for all intervals <= max_unpaired_length of targetRNA  *
* using the new RNAplfold of Vienna Package 1.8.x                             *
******************************************************************************/

void calc_ED_target_matrix_plfold(char* sequence)
{
   int mRNA_len = (int)strlen(sequence);
   float cutoff = 0.01;    /* bpcutoff for plfold */
   double ** pup = (double **) space((mRNA_len+1)*sizeof(double *)); //unpaired probs

   pup[0] = (double *)space(sizeof(double)); /*need only entry 0*/
   pup[0][0] = (double) max_unpaired_length;
   plist *dpp=NULL;

   pf_scale = -1;
   update_fold_params();
   pfl_fold(sequence, window_size, max_pair_dist, cutoff, pup, &dpp, NULL, NULL);

   // i+u: end position of unpaired region counting from 0, u: length of unpaired region
   for(int u=1; u<=max_unpaired_length; u++) {
      for(size_t i=0; i<mRNA_len-u+1; i++) {
         if (pup[i+u][u] == 0.0) {  // is not exact since we have double values, but pu can get really small
            ED_target[i][i+u-1] = 
               (fabs(ED_target_scaling - 0.0) < DOUBLE_DIFF) ? 0 : MAX_DOUBLE;
//            cout << "i,j: " << i <<"," << i+u-1 << " --> " << pup[i+u][u] << endl;
//            cout << "   " <<      " --> " << ED_target[i][i+u-1] << endl;
         }
         else {
            ED_target[i][i+u-1] = -RT * ED_target_scaling * log(pup[i+u][u]);
//            cout << "i,j: " << i <<"," << i+u-1 << " --> " << pup[i+u][u] << endl;
//            cout << "   " <<      " --> " << ED_target[i][i+u-1] << endl;
         }
      }
   }

    for (size_t i=0;i<=mRNA_len;i++) {
        free(pup[i]);
    }
    free(pup);
}


/******************************************************************************
* calculates ED values for all intervals <= max_unpaired_length of targetRNA  *
* using the old RNAplfold of Vienna Package 1.7.x                             *
******************************************************************************/

//void calc_ED_target_matrix_old_plfold(char* sequence)
//{
//   int len = (int)strlen(sequence);
////    double kT = (temperature+273.15)*1.98717/1000.0; /* in Kcal */
//   
//   float cutoff=0.01;
//   char* structure=NULL;
//   double *pup=NULL; /*prob of being unpaired*/
//   plist *pl;
//   
//   //cout << "RNAplfold used" << endl;
//      
//   for(int u=1;u<=max_unpaired_length;u++)
//   {
//      pf_scale = -1;
//      
//      pup=(double *)space((len+1)*sizeof(double));
//      pup[0]=(double)u;
//      
//      structure = (char *) space((unsigned) len+1);
//      update_fold_params();
//      
//      //printf("sequence: %s\n", sequence);
//      //printf("window_size: %d\n", window_size);
//      //printf("pairdist: %d\n", pairdist);
//      //printf("cutoff: %f\n", cutoff);
//      //cout << "pup: " << u << endl;
//      //cout << endl;
//      
//      pl=pfl_fold(sequence, window_size, max_pair_dist, cutoff, pup);
//      free(pl);
//      
//      pup_factor_correct(pup,len,u);
//      
//      //cout << endl << "u: " << u << endl;
//      for (int i=0; i<len-u+1; i++)
//      {
//         //bei pup nicht +1 sondern nun i+u, da das dort anscheinend so gespeichert wird (sieht zumindest im output so aus)
//         //(denn in pup[0] ist max_unpaired_length gespeichert)
//         if (pup[i+u] == 0.0)
//            ED_target[i][i+u-1] = MAX_DOUBLE;
//         else
//            ED_target[i][i+u-1] = -RT * log(pup[i+u]); 
//         //cout << "i: " << i << " --> " << pup[i+u] << endl;
//         //cout << "   " <<      " --> " << ED_target[i][i+u-1] << endl;
//         
//      }
//      //cout << endl;
//      free(structure);
//      free(pup);
//   }

//}


/*********************************************************************
* help function to scale pu_values calulated by the Vienna Package   *
* (this function is a modified version of a Vienna function)         *
*********************************************************************/

void pup_factor_correct(double *pup,int length,int unpaired_length) {

  float factor = 0;
  float tfact;
  bool left; //help variable to remind whether we are in left-case
  fflush(NULL);
  
  for (int i=unpaired_length; i<=length; i++) {
  
     left = false;
     //left
     if (i<=window_size) {
       factor=1./(i - unpaired_length + 1);
       left = true;
     }
     //right
     //if (i>length - window_size + unpaired_length - 1) {
     if (i>length - window_size + unpaired_length) {
     
       //if left == true, we are in left and right case at the same time and thus there exists no standard case
       //furthermore, tfact = 1/#of intervals of length winsize in length
       if (left == true)
          tfact = 1./(length - window_size + 1);
       else
          tfact=1./(length-i+1);
          
       if (tfact > factor) {
	  factor = tfact;
       }

     }
     //standard case
     else {
       tfact=1./(window_size - unpaired_length + 1);
       if (tfact > factor) {
	  factor = tfact;
       }
     }
     
     pup[i] *= factor;
     //printf("%d %.6g\n",i,pup[i]*factor);
  }


}


/******************************************************************************
* calculates PUvalue = exp(-(1/(kT))ED(start_unpaired,end_unfold))
******************************************************************************/

double calc_PU (double ED)
{
//    double kT = (temperature+273.15)*1.98717/1000.0; /* in Kcal */
   double PU = exp(-ED/RT);
   return PU;
}

double AU_penalty(char C1,char C2)
{
  if ((C1=='A' && C2=='U') ||
      (C1=='U' && C2=='A') ||
      (C1=='G' && C2=='U') ||
      (C1=='U' && C2=='G'))
    return terminal_AU;
  return 0.0;
}

// /*****************************************
// * 5' j   3'    5' di 3'
// * 3' id  5'    3'  j 5'
// *****************************************/
// double ed5(char i, char j, char d) 
// {
//   double energy=0;
//   if (i=='U' && j=='A')
//     {
//       if       (d=='A')   energy=-0.3;
//       else if  (d=='C')   energy=-0.1;
//       else if  (d=='G')   energy=-0.2;
//       else   /*(d=='U')*/ energy=-0.2;
//     }
//   else if (i=='G' && j=='C')
//     {
//       if       (d=='A')   energy=-0.2;
//       else if  (d=='C')   energy=-0.3;
//       else if  (d=='G')   energy=-0.0;
//       else   /*(d=='U')*/ energy=-0.0;
//     }
//   else if (i=='C' && j=='G')
//     {
//       if       (d=='A')   energy=-0.5;
//       else if  (d=='C')   energy=-0.3;
//       else if  (d=='G')   energy=-0.2;
//       else   /*(d=='U')*/ energy=-0.1;
//     }
//   else if (i=='U' && j=='G')
//     {
//       if       (d=='A')   energy=-0.3;
//       else if  (d=='C')   energy=-0.1;
//       else if  (d=='G')   energy=-0.2;
//       else   /*(d=='U')*/ energy=-0.2;
//     }
//   else if ( (i=='A' || i=='G') && j=='U')
//     {
//       if       (d=='A')   energy=-0.3;
//       else if  (d=='C')   energy=-0.3;
//       else if  (d=='G')   energy=-0.4;
//       else   /*(d=='U')*/ energy=-0.2;
//     }
//   return energy;
// }
// 
// /*****************************************
// * 5' jd 3'    5'  i 3'
// * 3' i  5'    3' dj 5'
// *****************************************/
// double ed3(char i, char j, char d) 
// {
//   double energy=0;
//   if (i=='U' && j=='A')
//     {
//       if       (d=='A')   energy=-0.8;
//       else if  (d=='C')   energy=-0.5;
//       else if  (d=='G')   energy=-0.8;
//       else   /*(d=='U')*/ energy=-0.6;
//     }
//   else if (i=='G' && j=='C')
//     {
//       if       (d=='A')   energy=-1.7;
//       else if  (d=='C')   energy=-0.8;
//       else if  (d=='G')   energy=-1.7;
//       else   /*(d=='U')*/ energy=-1.2;
//     }
//   else if (i=='C' && j=='G')
//     {
//       if       (d=='A')   energy=-1.1;
//       else if  (d=='C')   energy=-0.4;
//       else if  (d=='G')   energy=-1.3;
//       else   /*(d=='U')*/ energy=-0.6;
//     }
//   else if (i=='U' && j=='G')
//     {
//       if       (d=='A')   energy=-0.8;
//       else if  (d=='C')   energy=-0.5;
//       else if  (d=='G')   energy=-0.8;
//       else   /*(d=='U')*/ energy=-0.6;
//     }
//   else if ( (i=='A' || i=='G') && j=='U')
//     {
//       if       (d=='A')   energy=-0.7;
//       else if  (d=='C')   energy=-0.1;
//       else if  (d=='G')   energy=-0.7;
//       else   /*(d=='U')*/ energy=-0.1;
//     }
//   return energy;
// }


