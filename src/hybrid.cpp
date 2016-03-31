#include "hybrid.h"

/**********************************************************************************
*   allocates space for global matrices except ED_target and ED_miRNA             *
**********************************************************************************/

void allocate_global_matrices()
{
   // space allocation
   C                    = (double**) space(sizeof(double*)*strlen(sequence_target));
   Cseed                = (double**) space(sizeof(double*)*strlen(sequence_target));
   U                    = (double**) space(sizeof(double*)*strlen(sequence_target));
   che1                 = (int**)    space(sizeof(int*)   *strlen(sequence_target));
   che2                 = (int**)    space(sizeof(int*)   *strlen(sequence_target));
   ches1                = (int**)    space(sizeof(int*)   *strlen(sequence_target));
   ches2                = (int**)    space(sizeof(int*)   *strlen(sequence_target));
   seed                 = (double**) space(sizeof(double*)*strlen(sequence_target));
   seed_unpaired_target = (int**)    space(sizeof(int*)   *strlen(sequence_target));
   seed_unpaired_miRNA  = (int**)    space(sizeof(int*)   *strlen(sequence_target));
   hybrid_matrix        = (double*****) space(sizeof(double****)*strlen(sequence_target));

   for(int i=0;i<(int)strlen(sequence_target);i++)
   {
      C    [i]     = (double*) space(sizeof(double)*strlen(sequence_miRNA));
      Cseed[i]     = (double*) space(sizeof(double)*strlen(sequence_miRNA));
      U    [i]     = (double*) space(sizeof(double)*strlen(sequence_miRNA));
      che1 [i]     = (int*)    space(sizeof(int)   *strlen(sequence_miRNA));
      che2 [i]     = (int*)    space(sizeof(int)   *strlen(sequence_miRNA));
      ches1[i]     = (int*)    space(sizeof(int)   *strlen(sequence_miRNA));
      ches2[i]     = (int*)    space(sizeof(int)   *strlen(sequence_miRNA));
      seed [i]     = (double*) space(sizeof(double)*strlen(sequence_miRNA));
      seed_unpaired_target[i] = (int*) space(sizeof(int)*strlen(sequence_miRNA));
      seed_unpaired_miRNA [i] = (int*) space(sizeof(int)*strlen(sequence_miRNA));
      hybrid_matrix[i]        = (double****) space(sizeof(double***)*strlen(sequence_miRNA));
      for (int j=0; j<(int)strlen(sequence_miRNA); j++)
      {
         hybrid_matrix[i][j] = (double***) space(sizeof(double**)*(num_bp_seed+1));
         for (int k=0; k<=num_bp_seed; k++)
         {
            hybrid_matrix[i][j][k] = (double**) space(sizeof(double*)*(max_unpaired_target+1));
            for (int l=0; l<=max_unpaired_target; l++)
               hybrid_matrix[i][j][k][l] = (double*) space(sizeof(double)*(max_unpaired_miRNA+1));
         }
      }
   }
}

/**********************************************************************************
*   allocates space for ED_target matrix                                           *
**********************************************************************************/

void allocate_ED_target_matrix()
{
   // space allocation
   ED_target            = (double**) space(sizeof(double*)*strlen(sequence_target));

   for(int i=0;i<(int)strlen(sequence_target);i++)
   {
      ED_target[i] = (double*) space(sizeof(double)*strlen(sequence_target));
   }
}

/**********************************************************************************
*   allocates space for ED_miRNA matrix                                           *
**********************************************************************************/

void allocate_ED_miRNA_matrix()
{
   // space allocation
   ED_miRNA             = (double**) space(sizeof(double*)*strlen(sequence_miRNA ));

   for(int i=0;i<(int)strlen(sequence_miRNA);i++)
   {
      ED_miRNA[i] =(double*) space(sizeof(double)*strlen(sequence_miRNA));
   } 
}

/**********************************************************************************
*   initializing global matrices except ED_target and ED_miRNA                    *
**********************************************************************************/
  
void init_global_matrices()
{
  // init
  for(int i=0;i<(int)strlen(sequence_target);i++)
     for(int j=0;j<(int)strlen(sequence_miRNA);j++)
     {
	    C[i][j] = Cseed[i][j] = U[i][j] = seed[i][j] = MAX_DOUBLE;
	    che1[i][j] = che2[i][j] = ches1[i][j] = ches2[i][j] = -1;
	    seed_unpaired_target[i][j] = seed_unpaired_miRNA[i][j] = -1;
	    for (int k=0; k<=num_bp_seed; k++)
	       for (int l=0; l<=max_unpaired_target; l++)
	          for (int m=0; m<=max_unpaired_miRNA; m++)
	             hybrid_matrix[i][j][k][l][m] = MAX_DOUBLE;
     }    
}

/**********************************************************************************
*   initializing ED_target matrix                                                 *
**********************************************************************************/
  
void init_ED_target_matrix()
{
  // init
  for(int i=0;i<(int)strlen(sequence_target);i++)
    for(int j=0;j<(int)strlen(sequence_target);j++)
      ED_target[i][j] = MAX_DOUBLE;
}

/**********************************************************************************
*   initializing ED_miRNA matrix                                                 *
**********************************************************************************/
  
void init_ED_miRNA_matrix()
{
  // init
  for(int i=0;i<(int)strlen(sequence_miRNA);i++)
    for(int j=0;j<(int)strlen(sequence_miRNA);j++)
      ED_miRNA[i][j] = MAX_DOUBLE;
}

/**********************************************************************************
*       frees space for global matrices except ED_target and ED_miRNA             *
**********************************************************************************/

void free_global_matrices()
{
   for(int i=0;i<(int)strlen(sequence_target);i++)
   {
      free(C[i]);
      free(Cseed[i]);
      free(U[i]);
      free(che1[i]);
      free(che2[i]);
      free(ches1[i]);
      free(ches2[i]);
      free(seed[i]);
      free(seed_unpaired_target[i]);
      free(seed_unpaired_miRNA[i]);
	  
      for (int j=0; j<(int)strlen(sequence_miRNA); j++)
      {
         for (int k=0; k<=num_bp_seed; k++)
         {
            for (int l=0; l<=max_unpaired_target; l++)
               free(hybrid_matrix[i][j][k][l]);
			free(hybrid_matrix[i][j][k]);
         }
		 free(hybrid_matrix[i][j]);
      }
	  free(hybrid_matrix[i]);
   }
	  
   free(C);
   free(Cseed);
   free(U);
   free(che1);
   free(che2);
   free(ches1);
   free(ches2);
   free(seed);
   free(seed_unpaired_target);
   free(seed_unpaired_miRNA);
   free(hybrid_matrix);
}

/**********************************************************************************
*       frees space for ED_target matrix                                          *
**********************************************************************************/

void free_ED_target_matrix()
{
   for(int i=0;i<(int)strlen(sequence_target);i++)
   {
      free(ED_target[i]);
   }

   free(ED_target);
}

/**********************************************************************************
*       frees space for ED_miRNA matrix                                           *
**********************************************************************************/

void free_ED_miRNA_matrix()
{

   for(int i=0;i<(int)strlen(sequence_miRNA);i++)
   {
      free(ED_miRNA[i]);
   }

   free(ED_miRNA);
}

/**********************************************************************************
*      calculates one value of the hybrid_matrix                                  *
**********************************************************************************/

double calc_hybrid_value(int left_end_target, int left_end_miRNA, int num_bp, int unpaired_target, int unpaired_miRNA)
{
   double stacking, min_bulge_target, min_bulge_miRNA, min_interior_loop, bulge_target, bulge_miRNA, interior_loop;
   double min_energy = MAX_DOUBLE;

/*printf("left_end_target: %d\n", left_end_target);
printf("left_end_miRNA:  %d\n", left_end_miRNA);
printf("num_bp: %d\n", num_bp);
printf("unpaired_target: %d\n", unpaired_target);
printf("unpaired_miRNA:  %d\n", unpaired_miRNA);*/

   //test, whether both enclosing pairs can pair and enough bases are left on both sequences
   if ((left_end_target+num_bp+unpaired_target-1 >=(int)strlen(sequence_target)) ||
	(left_end_miRNA+num_bp+unpaired_miRNA-1 >= (int)strlen(sequence_miRNA)) ||
	(PairTest(sequence_target[left_end_target],sequence_miRNA[left_end_miRNA]) == false) ||
	(PairTest(sequence_target[left_end_target+num_bp+unpaired_target-1],sequence_miRNA[left_end_miRNA+num_bp+unpaired_miRNA-1]) == false))
   {
      return MAX_DOUBLE;
   }
   else
   {  //basic case
      if (num_bp == 2)
      {
         if ((unpaired_target == 0) && (unpaired_miRNA == 0)) /*stacking*/
            return StackingEnergy(int_sequence_target[left_end_target], int_sequence_miRNA[left_end_miRNA], int_sequence_target[left_end_target+1], int_sequence_miRNA[left_end_miRNA+1]);
            
         else if ((unpaired_target > 0) && (unpaired_target <= MAX_BULGE_SIZE) && (unpaired_miRNA == 0)) /*target_bulge, has to be at most MAX_BULGE_SIZE*/
            return BulgeEnergy(unpaired_target, int_sequence_target[left_end_target], int_sequence_miRNA[left_end_miRNA], int_sequence_target[left_end_target+1+unpaired_target], int_sequence_miRNA[left_end_miRNA+1]);
            
         else if ((unpaired_target == 0) && (unpaired_miRNA <= MAX_BULGE_SIZE) && (unpaired_miRNA > 0)) /*miRNA_bulge, has to be at most MAX_BULGE_SIZE*/
            return BulgeEnergy(unpaired_miRNA, int_sequence_target[left_end_target], int_sequence_miRNA[left_end_miRNA], int_sequence_target[left_end_target+1], int_sequence_miRNA[left_end_miRNA+1+unpaired_miRNA]);
            
         else if ((unpaired_target > 0) && (unpaired_target <= MAX_INTERIOR_ONE_SIDE_SIZE) && (unpaired_miRNA > 0) && (unpaired_miRNA <= MAX_INTERIOR_ONE_SIDE_SIZE)) /*interior loop, both sides have to be at most MAX_INTERIOR_ONE_SIDE_SIZE*/
            return InteriorLoopEnergy(unpaired_target, unpaired_miRNA, int_sequence_target[left_end_target], int_sequence_miRNA[left_end_miRNA], int_sequence_target[left_end_target+1+unpaired_target], int_sequence_miRNA[left_end_miRNA+1+unpaired_miRNA], int_sequence_target[left_end_target+1], int_sequence_miRNA[left_end_miRNA+1], int_sequence_target[left_end_target+unpaired_target], int_sequence_miRNA[left_end_miRNA+unpaired_miRNA]);   
         else
            return MAX_DOUBLE;
            
      }
      else if (num_bp > 2) //recursion
      {
         /*stacking*/
         stacking = Sum_or_MAX_DOUBLE(StackingEnergy(int_sequence_target[left_end_target], int_sequence_miRNA[left_end_miRNA], int_sequence_target[left_end_target+1], int_sequence_miRNA[left_end_miRNA+1]), hybrid_matrix[left_end_target+1][left_end_miRNA+1][num_bp-1][unpaired_target][unpaired_miRNA]);
         
         /* bulge in the target mRNA */         
         min_bulge_target = MAX_DOUBLE;
         for (int p=1; p<=Minimum(unpaired_target,MAX_BULGE_SIZE); p++)
         {
            bulge_target = Sum_or_MAX_DOUBLE(BulgeEnergy(p, int_sequence_target[left_end_target], int_sequence_miRNA[left_end_miRNA], int_sequence_target[left_end_target+1+p], int_sequence_miRNA[left_end_miRNA+1]), hybrid_matrix[left_end_target+1+p][left_end_miRNA+1][num_bp-1][unpaired_target-p][unpaired_miRNA]);
            
            if (bulge_target < min_bulge_target)
               min_bulge_target = bulge_target;
         }
         
         /* bulge in the miRNA */
         min_bulge_miRNA = MAX_DOUBLE;
         for (int q=1; q<=Minimum(unpaired_miRNA,MAX_BULGE_SIZE); q++)
         {
            bulge_miRNA = Sum_or_MAX_DOUBLE(BulgeEnergy(q, int_sequence_target[left_end_target], int_sequence_miRNA[left_end_miRNA], int_sequence_target[left_end_target+1], int_sequence_miRNA[left_end_miRNA+1+q]), hybrid_matrix[left_end_target+1][left_end_miRNA+1+q][num_bp-1][unpaired_target][unpaired_miRNA-q]);
            
            if (bulge_miRNA < min_bulge_miRNA)
               min_bulge_miRNA = bulge_miRNA;
         }
         
         /* interior loop */
         min_interior_loop = MAX_DOUBLE;
         for (int p=1; p<=Minimum(unpaired_target,MAX_INTERIOR_ONE_SIDE_SIZE); p++)
         {
            for (int q=1; q<=Minimum(unpaired_miRNA,MAX_INTERIOR_ONE_SIDE_SIZE); q++)
            {
               interior_loop = Sum_or_MAX_DOUBLE(InteriorLoopEnergy(p, q, int_sequence_target[left_end_target], int_sequence_miRNA[left_end_miRNA], int_sequence_target[left_end_target+1+p], int_sequence_miRNA[left_end_miRNA+1+q], int_sequence_target[left_end_target+1], int_sequence_miRNA[left_end_miRNA+1], int_sequence_target[left_end_target+p], int_sequence_miRNA[left_end_miRNA+q]), hybrid_matrix[left_end_target+1+p][left_end_miRNA+1+q][num_bp-1][unpaired_target-p][unpaired_miRNA-q]);
            
               if (interior_loop < min_interior_loop)
                  min_interior_loop = interior_loop;
            }
         }
         
         /* find the minimum of the 4 cases */
         min_energy = Minimum(min_energy, stacking);
         min_energy = Minimum(min_energy, min_bulge_target);
         min_energy = Minimum(min_energy, min_bulge_miRNA);
         min_energy = Minimum(min_energy, min_interior_loop);
         
         return min_energy;
      }
      else
      {
         cerr << "Forgotten hybrid cases!\n";
         return MAX_DOUBLE;
      }
   }
   return MAX_DOUBLE;
}



/********************************************************************************
* findes the best seed starting at base pair (i,j), having 'num_bp_seed' base   *
* pairs, having at most 'max_unpaired_target' unbound bases in the target, at   *
* most 'max_unpaired_miRNA' unbound bases in the miRNA, and at most             *
* 'max_unpaired' unbound bases at all                                           *
*                                                                               *
* the number of unbound bases in both sequences (which give the best seed) are  *
* stored in seed_unpaired_target(i,j) and seed_unpaired_miRNA(i,j)              *
********************************************************************************/

double compute_seed_ij(int i, int j)
{
   //test, whether it is a valid seed, i.e. j is within ncRNA seed region if applicable, x_i and y_j can pair and 
   //'num_bp_seed' base pairs starting at (i,j) are possible (depending on the lengths of the sequences)
   if ( (seed_region_ncRNA && (j<seed_region_start_ncRNA_3 || j+num_bp_seed-1>seed_region_end_ncRNA_3)) || 
		(PairTest(sequence_target[i],sequence_miRNA[j]) == false) ||
        (i+num_bp_seed-1 >= (int)strlen(sequence_target)) || (j+num_bp_seed-1 >= (int)strlen(sequence_miRNA)) )
   {
      return MAX_DOUBLE;
   }
   else
   {  

      double min_energy = MAX_DOUBLE;
      int min_target_unpaired, min_miRNA_unpaired; // numbers of unpaired bases when reaching the min_energy
      int i_star, j_star; //helper variables
      double hyb, energy; //energy helper during the recursion
      
      min_target_unpaired = min_miRNA_unpaired = -1;
      i_star = i+num_bp_seed-1;
      j_star = j+num_bp_seed-1;
      for (int unpaired_target = 0; unpaired_target <= max_unpaired_target; unpaired_target++)
         for (int unpaired_miRNA = 0; unpaired_miRNA <= max_unpaired_miRNA; unpaired_miRNA++)
         {
            if ((unpaired_target + unpaired_miRNA <= max_unpaired) && 
				(seed_region_ncRNA?(unpaired_miRNA <= seed_region_end_ncRNA_3-j_star):true))
            {

               hyb = hybrid_matrix[i][j][num_bp_seed][unpaired_target][unpaired_miRNA];

               // if first valid seed found, compute ED_target matrix, ED_miRNA matrix and C matrix if not done yet
               if (hyb < MAX_DOUBLE && (ED_miRNA_computed == false || ED_target_computed == false) ) {
                  if (ED_miRNA_computed == false) {
                    // calculate ED matrix for miRNA
                    if (use_RNAup) {
                       calc_ED_miRNA_matrix_up();
                     }
                    else {
                       calc_ED_miRNA_matrix();
                     }
                    ED_miRNA_computed = true;
                  }

                  if (ED_target_computed == false) {
                    // calculate ED matrix for target
                    if (use_RNAplfold)
                      calc_ED_target_matrix_plfold(sequence_target);
                    else
                    {
                      if (window)
                          calc_ED_target_matrix_window();
                      else
                          calc_ED_target_matrix(); 
                    }
                    ED_target_computed = true;
                 }

                  for(int i=(int)strlen(sequence_target)-1;i>=0;i--)
                     for(int j=(int)strlen(sequence_miRNA)-1;j>=0;j--)
                        C[i][j]=Cij(i,j);
               }

               if ((hyb < MAX_DOUBLE) && (C[i_star+unpaired_target][j_star+unpaired_miRNA] < MAX_DOUBLE) && 
					(ED_target[i][che1[i_star+unpaired_target][j_star+unpaired_miRNA]]<MAX_DOUBLE) && 
					(ED_miRNA[j][che2[i_star+unpaired_target][j_star+unpaired_miRNA]]<MAX_DOUBLE) 
					// next line is quick hack, please remove if not needed anymore, works only of unpaired_target and unpaired_miRNA is zero
					&& calc_PU(ED_target[i][i_star])*calc_PU(ED_miRNA[j][j_star]) >= pu_comb_threshold)
               {
                  energy = hyb + C[i_star+unpaired_target][j_star+unpaired_miRNA] - 
						   ED_target[i_star+unpaired_target][che1[i_star+unpaired_target][j_star+unpaired_miRNA]] -
						   ED_miRNA[j_star+unpaired_miRNA][che2[i_star+unpaired_target][j_star+unpaired_miRNA]] + 
						   ED_target[i][che1[i_star+unpaired_target][j_star+unpaired_miRNA]] + 
						   ED_miRNA[j][che2[i_star+unpaired_target][j_star+unpaired_miRNA]];
                  
                  if (energy < min_energy)
                  {
                     min_energy = energy;
                     min_target_unpaired = unpaired_target;
                     min_miRNA_unpaired = unpaired_miRNA;
                  }
               }
            }
         }

      seed_unpaired_target[i][j] = min_target_unpaired;
      seed_unpaired_miRNA[i][j] =  min_miRNA_unpaired;

      return min_energy;
   }
}





/****************************************************************************
* calculates the values of C[i][j]                                          *
****************************************************************************/

double Cij(int i,int j)
{
  /* test, whether bases can pair and are not the last bases in the sequences */
  if ( !PairTest(sequence_target[i],sequence_miRNA[j]) )
    return MAX_DOUBLE;
  else{
    double v1,v2,v3,v4,v5,v6,v7,v8;
    v1=v2=v3=v4=v5=v6=v7=v8=MAX_DOUBLE;
    
    /* stacking */
    v1=( (i+1>=(int)strlen(sequence_target) || j+1>=(int)strlen(sequence_miRNA) || C[i+1][j+1]==MAX_DOUBLE || ED_target[i][che1[i+1][j+1]]==MAX_DOUBLE || ED_miRNA[j][che2[i+1][j+1]]==MAX_DOUBLE)?MAX_DOUBLE:
	 StackingEnergy(int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1], int_sequence_miRNA[j+1])+
	 C[i+1][j+1]+
	 ED_target[i][che1[i+1][j+1]]-ED_target[i+1][che1[i+1][j+1]]+
	 ED_miRNA [j][che2[i+1][j+1]]-ED_miRNA [j+1][che2[i+1][j+1]]);
    
    /* bulge in the mRNA (target) */
    double min          = MAX_DOUBLE;
    int    best_bulge_k = 1;
    for(int k=1;k<=Minimum(MAX_BULGE_SIZE,(int)strlen(sequence_target)-i-2);k++)
      {
	min=( (i+1+k>=(int)strlen(sequence_target) || j+1>=(int)strlen(sequence_miRNA) || C[i+1+k][j+1]==MAX_DOUBLE || ED_target[i][che1[i+1+k][j+1]]==MAX_DOUBLE || ED_miRNA[j][che2[i+1+k][j+1]]==MAX_DOUBLE)?MAX_DOUBLE:
	      BulgeEnergy(k,int_sequence_target[i],int_sequence_miRNA[j],int_sequence_target[i+1+k],int_sequence_miRNA[j+1])+
	      C[i+1+k][j+1]+
	      ED_target[i][che1[i+1+k][j+1]]-ED_target[i+1+k][che1[i+1+k][j+1]]+
	      ED_miRNA [j][che2[i+1+k][j+1]]-ED_miRNA [j+1  ][che2[i+1+k][j+1]]);
	if (min<v2) 
	  { 
	    v2           = min;
	    best_bulge_k = k;
	  }
      }
    
    /* bulge in the miRNA */
    min              = MAX_DOUBLE;
    int best_bulge_l = 1;
    for(int l=1;l<=Minimum(MAX_BULGE_SIZE,(int)strlen(sequence_miRNA)-j-2);l++)
      {
	min=( (i+1>=(int)strlen(sequence_target) || j+1+l>=(int)strlen(sequence_miRNA) || C[i+1][j+1+l]==MAX_DOUBLE || ED_target[i][che1[i+1][j+1+l]]==MAX_DOUBLE || ED_miRNA[j][che2[i+1][j+1+l]]==MAX_DOUBLE)?MAX_DOUBLE:
	      BulgeEnergy(l,int_sequence_target[i],int_sequence_miRNA[j],int_sequence_target[i+1],int_sequence_miRNA[j+1+l])+
	      C[i+1][j+1+l]+
	      ED_target[i][che1[i+1][j+1+l]]-ED_target[i+1   ][che1[i+1][j+1+l]]+
	      ED_miRNA [j][che2[i+1][j+1+l]]-ED_miRNA [j+1+l ][che2[i+1][j+1+l]]);
	if (min<v3)
	  {
	    v3           = min;
	    best_bulge_l = l;
	  }
      }
    
    /* interior loop */
    min=MAX_DOUBLE;
    int best_il_k=1;
    int best_il_l=1;
    for(int k=1;k<=Minimum(MAX_INTERIOR_ONE_SIDE_SIZE,(int)strlen(sequence_target)-i-2);k++)
      for(int l=1;l<=Minimum(MAX_INTERIOR_ONE_SIDE_SIZE,(int)strlen(sequence_miRNA)-j-2);l++)
	{
	  min=( (i+1+k>=(int)strlen(sequence_target) || j+1+l>=(int)strlen(sequence_miRNA) || C[i+1+k][j+1+l]==MAX_DOUBLE || ED_target[i][che1[i+1+k][j+1+l]]==MAX_DOUBLE || ED_miRNA[j][che2[i+1+k][j+1+l]]==MAX_DOUBLE)?MAX_DOUBLE:
		InteriorLoopEnergy(k,l, int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1+k], int_sequence_miRNA[j+1+l], int_sequence_target[i+1], int_sequence_miRNA[j+1], int_sequence_target[i+k], int_sequence_miRNA[j+l])+ 
		C[i+1+k][j+1+l]+
		ED_target[i][che1[i+1+k][j+1+l]]-ED_target[i+1+k ][che1[i+1+k][j+1+l]]+
		ED_miRNA [j][che2[i+1+k][j+1+l]]-ED_miRNA [j+1+l ][che2[i+1+k][j+1+l]]);
	  if (min<v4)
	    {
	      v4        = min;
	      best_il_k = k;
	      best_il_l = l;
	    }
	}
    
    /* basic case: exterior loop */
//     v5 = ( (i+1>=(int)strlen(sequence_target) || j+1>=(int)strlen(sequence_miRNA)
//            || ED_target[i][i+1]==MAX_DOUBLE || ED_miRNA[j][j+1]==MAX_DOUBLE)?MAX_DOUBLE:
//            ed5(sequence_miRNA[j],sequence_target[i],sequence_miRNA[j+1])+
//            ed3(sequence_miRNA[j],sequence_target[i],sequence_target[i+1])+
//            ED_miRNA[j][j+1]+ED_target[i][i+1]+AU_penalty(sequence_miRNA[j],sequence_target[i]));
//     v6 = ( (i+1>=(int)strlen(sequence_target) || ED_target[i][i+1]==MAX_DOUBLE || ED_miRNA[j][j]==MAX_DOUBLE)?MAX_DOUBLE:
//            ed3(sequence_miRNA[j],sequence_target[i],sequence_target[i+1])+
//            ED_miRNA[j][j]+ED_target[i][i+1]+AU_penalty(sequence_miRNA[j],sequence_target[i]));
//     v7 = ( (j+1>=(int)strlen(sequence_miRNA)  || ED_target[i][i]==MAX_DOUBLE || ED_miRNA[j][j+1]==MAX_DOUBLE)?MAX_DOUBLE:
//            ed5(sequence_miRNA[j],sequence_target[i],sequence_miRNA[j+1])+
//            ED_miRNA[j][j+1]+ED_target[i][i]+AU_penalty(sequence_miRNA[j],sequence_target[i]));
//     v8 = ( (ED_target[i][i] == MAX_DOUBLE || ED_miRNA[j][j] == MAX_DOUBLE)?MAX_DOUBLE:
//            (AU_penalty(sequence_miRNA[j],sequence_target[i])+ED_target[i][i]+ED_miRNA[j][j]));

    v5 = ( (i+1>=(int)strlen(sequence_target) || j+1>=(int)strlen(sequence_miRNA)
           || ED_target[i][i+1]==MAX_DOUBLE || ED_miRNA[j][j+1]==MAX_DOUBLE)?MAX_DOUBLE:
           ed5(int_sequence_miRNA[j],int_sequence_target[i],int_sequence_miRNA[j+1])+
           ed3(int_sequence_miRNA[j],int_sequence_target[i],int_sequence_target[i+1])+
           ED_miRNA[j][j+1]+ED_target[i][i+1]+AU_penalty(sequence_miRNA[j],sequence_target[i])+duplex_init);
    v6 = ( (i+1>=(int)strlen(sequence_target) || ED_target[i][i+1]==MAX_DOUBLE || ED_miRNA[j][j]==MAX_DOUBLE)?MAX_DOUBLE:
           ed3(int_sequence_miRNA[j],int_sequence_target[i],int_sequence_target[i+1])+
           ED_miRNA[j][j]+ED_target[i][i+1]+AU_penalty(sequence_miRNA[j],sequence_target[i])+duplex_init);
    v7 = ( (j+1>=(int)strlen(sequence_miRNA)  || ED_target[i][i]==MAX_DOUBLE || ED_miRNA[j][j+1]==MAX_DOUBLE)?MAX_DOUBLE:
           ed5(int_sequence_miRNA[j],int_sequence_target[i],int_sequence_miRNA[j+1])+
           ED_miRNA[j][j+1]+ED_target[i][i]+AU_penalty(sequence_miRNA[j],sequence_target[i])+duplex_init);
    v8 = ( (ED_target[i][i] == MAX_DOUBLE || ED_miRNA[j][j] == MAX_DOUBLE)?MAX_DOUBLE:
           (AU_penalty(sequence_miRNA[j],sequence_target[i])+ED_target[i][i]+ED_miRNA[j][j])+duplex_init);

    // che values
   if (heuristic) {
      if      (v1<=v2 && v1<=v3 && v1<=v4 && v1<=v5 && v1<=v6 && v1<=v7 && v1<=v8 && v1<MAX_DOUBLE) { 
                  che1[i][j]=che1[i+1][j+1]; che2[i][j]=che2[i+1][j+1]; return v1 ;} 
      else if (v2<=v1 && v2<=v3 && v2<=v4 && v2<=v5 && v2<=v6 && v2<=v7 && v2<=v8 && v2<MAX_DOUBLE) {
                  che1[i][j]=che1[i+1+best_bulge_k][j+1]; che2[i][j]=che2[i+1+best_bulge_k][j+1]; return v2 ;} 
      else if (v3<=v1 && v3<=v2 && v3<=v4 && v3<=v5 && v3<=v6 && v3<=v7 && v3<=v8 && v3<MAX_DOUBLE) {
                  che1[i][j]=che1[i+1][j+1+best_bulge_l]; che2[i][j]=che2[i+1][j+1+best_bulge_l]; return v3 ;} 
      else if (v4<=v1 && v4<=v2 && v4<=v3 && v4<=v5 && v4<=v6 && v4<=v7 && v4<=v8 && v4<MAX_DOUBLE) {
                  che1[i][j]=che1[i+1+best_il_k][j+1+best_il_l]; che2[i][j]=che2[i+1+best_il_k][j+1+best_il_l]; return v4 ;}
      else if (v5<=v1 && v5<=v2 && v5<=v3 && v5<=v4 && v5<=v6 && v5<=v7 && v5<=v8 && v5<MAX_DOUBLE) {
                  che1[i][j]=i+1; che2[i][j]=j+1; return v5 ;}
      else if (v6<=v1 && v6<=v2 && v6<=v3 && v6<=v4 && v6<=v5 && v6<=v7 && v6<=v8 && v6<MAX_DOUBLE) {
                  che1[i][j]=i+1; che2[i][j]=j; return v6 ;}
      else if (v7<=v1 && v7<=v2 && v7<=v3 && v7<=v4 && v7<=v5 && v7<=v6 && v7<=v8 && v7<MAX_DOUBLE) {
                  che1[i][j]=i; che2[i][j]=j+1; return v7 ;}
      else if (v8<MAX_DOUBLE) { che1[i][j]=i; che2[i][j]=j; return v8 ;}
   }
   // without heuristic keep che values fixed
   else {
      // allow end of hybridization only if che value is set to current value of i and j
      if (v5<=v1 && v5<=v2 && v5<=v3 && v5<=v4 && v5<=v6 && v5<=v7 && v5<=v8 && v5<MAX_DOUBLE
               && i+1 == che1[i][j] && j+1 == che2[i][j]) {
                  return v5 ;}
      else if (v6<=v1 && v6<=v2 && v6<=v3 && v6<=v4 && v6<=v5 && v6<=v7 && v6<=v8 && v6<MAX_DOUBLE
               && i+1 == che1[i][j] && j == che2[i][j]) {
                  return v6 ;}
      else if (v7<=v1 && v7<=v2 && v7<=v3 && v7<=v4 && v7<=v5 && v7<=v6 && v7<=v8 && v7<MAX_DOUBLE
               && i == che1[i][j] && j+1 == che2[i][j]) {
                  return v7 ;}
      else if (v8<=v1 && v8<=v2 && v8<=v3 && v8<=v4 && v8<=v5 && v8<=v6 && v8<=v7 && v8<MAX_DOUBLE
               && i == che1[i][j] && j == che2[i][j]) { 
                  return v8 ;}
      // v5 to v8 are not valid if che value is not set to current value of i and j
      else if (v1<=v2 && v1<=v3 && v1<=v4 && v1<MAX_DOUBLE) { return v1 ;} 
      else if (v2<=v1 && v2<=v3 && v2<=v4 && v2<MAX_DOUBLE) { return v2 ;} 
      else if (v3<=v1 && v3<=v2 && v3<=v4 && v3<MAX_DOUBLE) { return v3 ;} 
      else if (v4<=v1 && v4<=v2 && v4<=v3 && v4<MAX_DOUBLE) { return v4 ;}

   }
  }
  return MAX_DOUBLE;
}


/****************************************************************************
* calculates the values of C^seed[i][j]                                     *
****************************************************************************/

double Cij_seed(int i,int j)
{
  /* test, whether bases can pair */
  if ( !PairTest(sequence_target[i],sequence_miRNA[j]) )
    return MAX_DOUBLE;
  else{
    double v1,v2,v3,v4,v5;
    v1=v2=v3=v4=v5=MAX_DOUBLE;

    /* stacking */
    v1=( (i+1>=(int)strlen(sequence_target) || j+1>=(int)strlen(sequence_miRNA) || Cseed[i+1][j+1]==MAX_DOUBLE || ED_target[i][ches1[i+1][j+1]]==MAX_DOUBLE || ED_miRNA[j][ches2[i+1][j+1]]==MAX_DOUBLE)?MAX_DOUBLE:
	 StackingEnergy(int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1], int_sequence_miRNA[j+1])+
	 Cseed[i+1][j+1]+
	 ED_target[i][ches1[i+1][j+1]]-ED_target[i+1][ches1[i+1][j+1]]+
	 ED_miRNA [j][ches2[i+1][j+1]]-ED_miRNA [j+1][ches2[i+1][j+1]]);

    /* bulge in the mRNA (target) */
    double min          = MAX_DOUBLE;
    int    best_bulge_k = 1;
    for(int k=1;k<=Minimum(MAX_BULGE_SIZE,(int)strlen(sequence_target)-i-2);k++)
      {
	min=( (i+1+k>=(int)strlen(sequence_target) || j+1>=(int)strlen(sequence_miRNA) || Cseed[i+1+k][j+1]==MAX_DOUBLE || ED_target[i][ches1[i+1+k][j+1]]==MAX_DOUBLE || ED_miRNA[j][ches2[i+1+k][j+1]]==MAX_DOUBLE)?MAX_DOUBLE:
	     BulgeEnergy(k,int_sequence_target[i],int_sequence_miRNA[j],int_sequence_target[i+1+k],int_sequence_miRNA[j+1])+
	     Cseed[i+1+k][j+1]+
	     ED_target[i][ches1[i+1+k][j+1]]-ED_target[i+1+k][ches1[i+1+k][j+1]]+
	     ED_miRNA [j][ches2[i+1+k][j+1]]-ED_miRNA [j+1  ][ches2[i+1+k][j+1]]);
	if (min<v2) 
	  { 
	    v2           = min;
	    best_bulge_k = k;
	  }
      }

    /* bulge in the miRNA */
    min              = MAX_DOUBLE;
    int best_bulge_l = 1;
    for(int l=1;l<=Minimum(MAX_BULGE_SIZE,(int)strlen(sequence_miRNA)-j-2);l++)
      {
	min=( (i+1>=(int)strlen(sequence_target) || j+1+l>=(int)strlen(sequence_miRNA) || Cseed[i+1][j+1+l]==MAX_DOUBLE || ED_target[i][ches1[i+1][j+1+l]]==MAX_DOUBLE || ED_miRNA[j][ches2[i+1][j+1+l]]==MAX_DOUBLE)?MAX_DOUBLE:
	     BulgeEnergy(l,int_sequence_target[i],int_sequence_miRNA[j],int_sequence_target[i+1],int_sequence_miRNA[j+1+l])+
	     Cseed[i+1][j+1+l]+
	     ED_target[i][ches1[i+1][j+1+l]]-ED_target[i+1   ][ches1[i+1][j+1+l]]+
	     ED_miRNA [j][ches2[i+1][j+1+l]]-ED_miRNA [j+1+l ][ches2[i+1][j+1+l]]);
	if (min<v3)
	  {
	    v3           = min;
	    best_bulge_l = l;
	  }
      }

    /* interior loop */
    min           = MAX_DOUBLE;
    int best_il_k = 1;
    int best_il_l = 1;
    for(int k=1;k<=Minimum(MAX_INTERIOR_ONE_SIDE_SIZE,(int)strlen(sequence_target)-i-2);k++)
      for(int l=1;l<=Minimum(MAX_INTERIOR_ONE_SIDE_SIZE,(int)strlen(sequence_miRNA)-j-2);l++)
	{
	  min=( (i+1+k>=(int)strlen(sequence_target) || j+1+l>=(int)strlen(sequence_miRNA) || Cseed[i+1+k][j+1+l]==MAX_DOUBLE || ED_target[i][ches1[i+1+k][j+1+l]]==MAX_DOUBLE || ED_miRNA[j][ches2[i+1+k][j+1+l]]==MAX_DOUBLE)?MAX_DOUBLE:
	       InteriorLoopEnergy(k,l, int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1+k], int_sequence_miRNA[j+1+l], int_sequence_target[i+1], int_sequence_miRNA[j+1], int_sequence_target[i+k], int_sequence_miRNA[j+l])+ 
	       Cseed[i+1+k][j+1+l]+
	       ED_target[i][ches1[i+1+k][j+1+l]]-ED_target[i+1+k ][ches1[i+1+k][j+1+l]]+
	       ED_miRNA [j][ches2[i+1+k][j+1+l]]-ED_miRNA [j+1+l ][ches2[i+1+k][j+1+l]]);
	  if (min<v4)
	    {
	      v4        = min;
	      best_il_k = k;
	      best_il_l = l;
	    }
	}
    
    /* basic case: exterior loop */
    v5=seed[i][j];
    
    // ches values
    if (v1<=v2 && v1<=v3 && v1<=v4 && v1<=v5 && v1<MAX_DOUBLE) 
    {  
       ches1[i][j]=ches1[i+1][j+1]; 
       ches2[i][j]=ches2[i+1][j+1]; 
       return v1 ;
    } 
    else if (v2<=v1 && v2<=v3 && v2<=v4 && v2<=v5 && v2<MAX_DOUBLE) 
    { 
       ches1[i][j]=ches1[i+1+best_bulge_k][j+1]; 
       ches2[i][j]=ches2[i+1+best_bulge_k][j+1]; 
       return v2 ;
    } 
    else if (v3<=v1 && v3<=v2 && v3<=v4 && v3<=v5 && v3<MAX_DOUBLE) 
    { 
       ches1[i][j]=ches1[i+1][j+1+best_bulge_l]; 
       ches2[i][j]=ches2[i+1][j+1+best_bulge_l]; 
       return v3 ;
    } 
    else if (v4<=v1 && v4<=v2 && v4<=v3 && v4<=v5 && v4<MAX_DOUBLE) 
    { 
       ches1[i][j]=ches1[i+1+best_il_k][j+1+best_il_l];
       ches2[i][j]=ches2[i+1+best_il_k][j+1+best_il_l]; 
       return v4 ;
    }
    else if (v5 < MAX_DOUBLE) /* seed found at pos (i,j)*/
    {
       ches1[i][j] = che1[i+num_bp_seed+seed_unpaired_target[i][j]-1][j+num_bp_seed+seed_unpaired_miRNA[i][j]-1]; 
       ches2[i][j] = che2[i+num_bp_seed+seed_unpaired_target[i][j]-1][j+num_bp_seed+seed_unpaired_miRNA[i][j]-1];
       return v5 ;
    }
  }
  /* no seed found (v5 == MAX_DOUBLE) and all other values (v1, v2, v3, v4 == MAX_DOUBLE) */
  return MAX_DOUBLE;
}

double U_ij(int i,int j)
{

  if (i<(int)strlen(sequence_target)-1 && j<(int)strlen(sequence_miRNA)-1) {
      /*
      cout << "*************************" << endl;
      cout << "  i : " << i << " j : " << j << endl;
      cout << "U_ij 1 : " << U[i+1][j] << endl;
      cout << "U_ij 2 : " << U[i][j+1] << endl;
      cout << "U_ij 3 : " << (Cseed[i+1][j+1]==MAX_DOUBLE?MAX_DOUBLE:Cseed[i+1][j+1]+ed5(sequence_target[i+1],sequence_miRNA[j+1],sequence_target[i])+ed3(sequence_target[i+1],sequence_miRNA[j+1],sequence_miRNA [j])+AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])) << " " 
	                  << Cseed[i+1][j+1] << " " 
	                  <<(Cseed[i+1][j+1]==MAX_DOUBLE) << " " << endl;
      */
    return Minimum(
            U[i+1][j],Minimum(
               U[i][j+1],Cseed[i+1][j+1]==MAX_DOUBLE?MAX_DOUBLE:Minimum(
                  Cseed[i+1][j+1]+ed5(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_target[i])+ed3(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_miRNA[j])
                  +AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])+ED_target[i][ches1[i+1][j+1]]+ED_miRNA[j][ches2[i+1][j+1]]
                  -ED_target[i+1][ches1[i+1][j+1]]-ED_miRNA[j+1][ches2[i+1][j+1]],Minimum(
                     Cseed[i+1][j+1]+ed5(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_target[i])+AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])
                     +ED_target[i][ches1[i+1][j+1]]-ED_target[i+1][ches1[i+1][j+1]],Minimum(
                        Cseed[i+1][j+1]+ed3(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_miRNA[j])+AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])
                        +ED_miRNA[j][ches2[i+1][j+1]]-ED_miRNA[j+1][ches2[i+1][j+1]],
                        Cseed[i+1][j+1]+AU_penalty(sequence_target[i+1],sequence_miRNA[j+1]))))));
    }
  // cases U(miRNA_length,target_length), U(miRNA_length-1,target_length)
  // and U(miRNA_length,target_length-1) which means no hybridization
  else
    return 0.0;
}


/****************************************************************************
* finds the best hybridizing part                                           *
****************************************************************************/

double compute_hybrid()
{
  double v, min1, min2, hybrid_energy;
  v = min1 = min2 = hybrid_energy = MAX_DOUBLE;

  // compute recurrences
  if (heuristic) {
  //has been moved to compute_seed_ij and is computed only when valid seed has been found
//   for(int i=(int)strlen(sequence_target)-1;i>=0;i--)
//     for(int j=(int)strlen(sequence_miRNA)-1;j>=0;j--)
//       C[i][j]=Cij(i,j);

  for (int i=(int)strlen(sequence_target)-1;i>=0;i--)
    for(int j=(int)strlen(sequence_miRNA)-1;j>=0;j--)
       for (int l=2; l<=num_bp_seed; l++)
          for (int b_target=0; b_target<=max_unpaired_target; b_target++)
             for (int b_miRNA=0; b_miRNA<=max_unpaired_miRNA; b_miRNA++)
                 hybrid_matrix[i][j][l][b_target][b_miRNA] = calc_hybrid_value(i,j,l,b_target,b_miRNA);

  for(int i=(int)strlen(sequence_target)-1;i>=0;i--)
    for(int j=(int)strlen(sequence_miRNA)-1;j>=0;j--)
      {
	seed [i][j] = compute_seed_ij(i,j);
        //is now computed only when valid seed has been found
// 	Cseed[i][j] = Cij_seed(i,j);
// 	U    [i][j] = U_ij(i,j);
      }

  //is now computed only when valid seed has been found
  if (ED_target_computed && ED_miRNA_computed) {
    for(int i=(int)strlen(sequence_target)-1;i>=0;i--)
        for(int j=(int)strlen(sequence_miRNA)-1;j>=0;j--)
        {
            Cseed[i][j] = Cij_seed(i,j);
            U    [i][j] = U_ij(i,j);
        }
  }

  for (int i=1; i<(int)strlen(sequence_target);i++)
    {
      v = ( Cseed[i][0]==MAX_DOUBLE ? MAX_DOUBLE : Minimum(Cseed[i][0]+AU_penalty(sequence_target[i],sequence_miRNA[0]),
             Cseed[i][0]+ed5(int_sequence_target[i],int_sequence_miRNA[0],int_sequence_target[i-1])+AU_penalty(sequence_target[i],sequence_miRNA[0])
             +ED_target[i-1][ches1[i][0]]-ED_target[i][ches1[i][0]] ) );
      if (v < min1)
        {
          min1 = v;
        }
    }

  for (int j=1; j<(int)strlen(sequence_miRNA);j++)
    {
      v = ( Cseed[0][j]==MAX_DOUBLE ? MAX_DOUBLE : Minimum(Cseed[0][j]+AU_penalty(sequence_target[0],sequence_miRNA[j]),
             Cseed[0][j]+ed3(int_sequence_target[0],int_sequence_miRNA[j],int_sequence_miRNA[j-1])+AU_penalty(sequence_target[0],sequence_miRNA[j])
             +ED_miRNA[j-1][ches2[0][j]]-ED_miRNA[j][ches2[0][j]] ) );
      if (v < min2)
        {
          min2 = v;
        }
    }

  hybrid_energy= Minimum(0.0,
                  Minimum(Cseed[0][0]==MAX_DOUBLE?MAX_DOUBLE:Cseed[0][0]+AU_penalty(sequence_target[0],sequence_miRNA[0]),
                     Minimum(min1,
                        Minimum(min2, U[0][0]))));
  }

  else {
  double hybrid_energy_pq = MAX_DOUBLE;

   // calculate ED matrix for miRNA
   if (use_RNAup)
      calc_ED_miRNA_matrix_up();
   else
      calc_ED_miRNA_matrix();
   ED_miRNA_computed = true;

   // calculate ED matrix for target
   if (use_RNAplfold)
       calc_ED_target_matrix_plfold(sequence_target);
   else {
      if (window)
         calc_ED_target_matrix_window();
      else
         calc_ED_target_matrix(); 
   }
   ED_target_computed = true;

  // compute recurrences without heuristic for each possible hybridization end p,q
  for(int p=(int)strlen(sequence_target)-1;p>=0;p--)
    for(int q=(int)strlen(sequence_miRNA)-1;q>=0;q--) {
      // set che values to p,q
      init_global_matrices();
      for(int i=0;i<=p;i++)
          for(int j=0;j<=q;j++) {
            che1[i][j] = p;
            che2[i][j] = q;
          }

        // compute matrices for all indices <= p,q
        for(int i=p;i>=0;i--)
           for(int j=q;j>=0;j--)
             C[i][j]=Cij(i,j);

        for (int i=p;i>=0;i--)
          for(int j=q;j>=0;j--)
            for (int l=2; l<=num_bp_seed; l++)
              for (int b_target=0; b_target<=max_unpaired_target; b_target++)
                for (int b_miRNA=0; b_miRNA<=max_unpaired_miRNA; b_miRNA++)
                  hybrid_matrix[i][j][l][b_target][b_miRNA] = calc_hybrid_value(i,j,l,b_target,b_miRNA);

        for(int i=p;i>=0;i--) {
          for(int j=q;j>=0;j--) {
            seed [i][j] = compute_seed_ij(i,j);
            Cseed[i][j] = Cij_seed(i,j);
            U    [i][j] = U_ij(i,j);
          }
        }

        for (int i=1; i<(int)strlen(sequence_target);i++) {
          v = ( Cseed[i][0]==MAX_DOUBLE ? MAX_DOUBLE : Minimum(Cseed[i][0]+AU_penalty(sequence_target[i],sequence_miRNA[0]),
                 Cseed[i][0]+ed5(int_sequence_target[i],int_sequence_miRNA[0],int_sequence_target[i-1])+AU_penalty(sequence_target[i],sequence_miRNA[0])
                 +ED_target[i-1][ches1[i][0]]-ED_target[i][ches1[i][0]] ) );
          if (v < min1) {
            min1 = v;
          }
        }
        for (int j=1; j<(int)strlen(sequence_miRNA);j++) {
          v = ( Cseed[0][j]==MAX_DOUBLE ? MAX_DOUBLE : Minimum(Cseed[0][j]+AU_penalty(sequence_target[0],sequence_miRNA[j]),
                 Cseed[0][j]+ed3(int_sequence_target[0],int_sequence_miRNA[j],int_sequence_miRNA[j-1])+AU_penalty(sequence_target[0],sequence_miRNA[j])
                 +ED_miRNA[j-1][ches2[0][j]]-ED_miRNA[j][ches2[0][j]] ) );
          if (v < min2) {
            min2 = v;
          }
        }

        hybrid_energy_pq = Minimum(0.0,
                              Minimum(Cseed[0][0]==MAX_DOUBLE?MAX_DOUBLE:Cseed[0][0]+AU_penalty(sequence_target[0],sequence_miRNA[0]),
                                 Minimum(min1,
                                    Minimum(min2, U[0][0]))));

        if (hybrid_energy_pq < hybrid_energy) {
          best_che1 = p;
          best_che2 = q;
          hybrid_energy = hybrid_energy_pq;
        }
      }
  }

  /*
  cout << "res 1 : " << 0 << endl;
  cout << "res 2 : " << (Cseed[0][0]==MAX_DOUBLE?MAX_DOUBLE:Cseed[0][0]+AU_penalty(sequence_target[0],sequence_miRNA[0])) << endl;
  cout << "res 3 : " << (Cseed[1][0]==MAX_DOUBLE?MAX_DOUBLE:Cseed[1][0]+ed5(sequence_target[1],sequence_miRNA[0],sequence_target[0])+AU_penalty(sequence_target[1],sequence_miRNA[0])) << endl;
  cout << "res 4 : " << (Cseed[0][1]==MAX_DOUBLE?MAX_DOUBLE:Cseed[0][1]+ed3(sequence_target[0],sequence_miRNA[1],sequence_miRNA[0])+AU_penalty(sequence_target[0],sequence_miRNA[1]))  << endl;
  cout << "res 5 : " << U[0][0] << endl;
  */

   return hybrid_energy;
}

/****************************************************************************
* finds best hybridization at positions not constrained by                  *
* hybrid_pos_target                                                         *
* calculates traceback of suboptimal hybridization                          *
* returns whether such a suboptimal hybridization was found and the energy  *
* of that hybridization (if applicable)                                     *
****************************************************************************/

pair<bool, double>  find_next_subopt() {
  bool hybridization_found = false;
  double hybrid_energy = MAX_DOUBLE;
  // hybridization energy of 0.0 means no significant hybridization
  double min_hybrid_energy = 0.0;
  int hybrid_start_target, hybrid_start_miRNA;
  hybrid_start_target = hybrid_start_miRNA = -1;

  // test position (0,0) for being a hybridization not considered yet
  if (ches1[0][0] != -1 && hybrid_pos_target.substr(0,ches1[0][0]+1).find('1',0) == string::npos &&
       (Cseed[0][0]==MAX_DOUBLE ? MAX_DOUBLE :
       (hybrid_energy=Cseed[0][0]+AU_penalty(sequence_target[0],sequence_miRNA[0]))) < min_hybrid_energy)
   {
     min_hybrid_energy = hybrid_energy;
     hybridization_found = true;
     hybrid_start_target = 0;
     hybrid_start_miRNA = 0;
     pos_target_1_dangle=hybrid_start_target;
     pos_miRNA_1_dangle =hybrid_start_miRNA;
   }

  // test position (i,0) for being a hybridization not considered yet
  for(int i=1; i<(int)strlen(sequence_target);i++) {
     if (ches1[i][0] != -1 && hybrid_pos_target.substr(i,ches1[i][0]-i+1).find('1',0) == string::npos &&
          (Cseed[i][0]==MAX_DOUBLE ? MAX_DOUBLE : 
          (hybrid_energy=Cseed[i][0]+AU_penalty(sequence_target[i],sequence_miRNA[0]))) < min_hybrid_energy)
      {
       min_hybrid_energy = hybrid_energy;
       hybridization_found = true;
       hybrid_start_target = i;
       hybrid_start_miRNA = 0;
       pos_target_1_dangle=hybrid_start_target;
       pos_miRNA_1_dangle =hybrid_start_miRNA;
      }
     if (ches1[i][0] != -1 && hybrid_pos_target.substr(i,ches1[i][0]-i+1).find('1',0) == string::npos &&
          (Cseed[i][0]==MAX_DOUBLE ? MAX_DOUBLE : 
          (hybrid_energy=Cseed[i][0]+ed5(int_sequence_target[i],int_sequence_miRNA[0],int_sequence_target[i-1])
						 +AU_penalty(sequence_target[i],sequence_miRNA[0])+ED_target[i-1][ches1[i][0]]-ED_target[i][ches1[i][0]]))
						 < min_hybrid_energy)
      {
       min_hybrid_energy = hybrid_energy;
       hybridization_found = true;
       hybrid_start_target = i;
       hybrid_start_miRNA = 0;
       pos_target_1_dangle=hybrid_start_target-1;
       pos_miRNA_1_dangle =hybrid_start_miRNA;
      }
  }

  // test position (0,j) for being a hybridization not considered yet
  for(int j=1; j<(int)strlen(sequence_miRNA);j++) {
     if (ches1[0][j] != -1 && hybrid_pos_target.substr(0,ches1[0][j]+1).find('1',0) == string::npos &&
          (Cseed[0][j]==MAX_DOUBLE ? MAX_DOUBLE :
          (hybrid_energy=Cseed[0][j]+AU_penalty(sequence_target[0],sequence_miRNA[j]))) < min_hybrid_energy)
      {
       min_hybrid_energy = hybrid_energy;
       hybridization_found = true;
       hybrid_start_target = 0;
       hybrid_start_miRNA = j;
       pos_target_1_dangle=hybrid_start_target;
       pos_miRNA_1_dangle =hybrid_start_miRNA;
      }
     if (ches1[0][j] != -1 && hybrid_pos_target.substr(0,ches1[0][j]+1).find('1',0) == string::npos &&
          (Cseed[0][j]==MAX_DOUBLE ? MAX_DOUBLE :
          (hybrid_energy=Cseed[0][j]+ed3(int_sequence_target[0],int_sequence_miRNA[j],int_sequence_miRNA[j-1])
						 +AU_penalty(sequence_target[0],sequence_miRNA[j])+ED_miRNA[j-1][ches2[0][j]]-ED_miRNA[j][ches2[0][j]]))
						 < min_hybrid_energy)
      {
       min_hybrid_energy = hybrid_energy;
       hybridization_found = true;
       hybrid_start_target = 0;
       hybrid_start_miRNA = j;
       pos_target_1_dangle=hybrid_start_target;
       pos_miRNA_1_dangle =hybrid_start_miRNA-1;
      }
  }

  // test positions (i,j) for being a hybridization not considered yet
  for(int i=0; i<(int)strlen(sequence_target)-1;i++) {
    for(int j=0; j<(int)strlen(sequence_miRNA)-1;j++) {
       if (ches1[i+1][j+1] != -1 && hybrid_pos_target.substr(i+1,ches1[i+1][j+1]-i).find('1',0) == string::npos &&
           (Cseed[i+1][j+1]==MAX_DOUBLE ? MAX_DOUBLE : 
           (hybrid_energy=Cseed[i+1][j+1]+ed5(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_target[i])
                          +ed3(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_miRNA[j])+AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])
                          +ED_target[i][ches1[i+1][j+1]]+ED_miRNA[j][ches2[i+1][j+1]]-ED_target[i+1][ches1[i+1][j+1]]-ED_miRNA[j+1][ches2[i+1][j+1]]))
			< min_hybrid_energy)
        {
         min_hybrid_energy = hybrid_energy;
         hybridization_found = true;
         hybrid_start_target = i+1;
         hybrid_start_miRNA = j+1;
		 pos_target_1_dangle = i;
         pos_miRNA_1_dangle =j;
		}
		if (ches1[i+1][j+1] != -1 && hybrid_pos_target.substr(i+1,ches1[i+1][j+1]-i).find('1',0) == string::npos &&
           (Cseed[i+1][j+1]==MAX_DOUBLE ? MAX_DOUBLE : 
           (hybrid_energy=Cseed[i+1][j+1]+ed5(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_target[i])
						  +AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])+ED_target[i][ches1[i+1][j+1]]-ED_target[i+1][ches1[i+1][j+1]]))
			< min_hybrid_energy)
        {
         min_hybrid_energy = hybrid_energy;
         hybridization_found = true;
         hybrid_start_target = i+1;
         hybrid_start_miRNA = j+1;
         pos_target_1_dangle = i;
         pos_miRNA_1_dangle =j+1;
		}
		if (ches1[i+1][j+1] != -1 && hybrid_pos_target.substr(i+1,ches1[i+1][j+1]-i).find('1',0) == string::npos &&
           (Cseed[i+1][j+1]==MAX_DOUBLE ? MAX_DOUBLE : 
           (hybrid_energy=Cseed[i+1][j+1]+ed3(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_miRNA[j])
						  +AU_penalty(sequence_target[i+1],sequence_miRNA[j+1]) +ED_miRNA[j][ches2[i+1][j+1]]-ED_miRNA[j+1][ches2[i+1][j+1]]))
			< min_hybrid_energy)
        {
         min_hybrid_energy = hybrid_energy;
         hybridization_found = true;
         hybrid_start_target = i+1;
         hybrid_start_miRNA = j+1;
         pos_target_1_dangle = i+1;
         pos_miRNA_1_dangle =j;
		}
		if (ches1[i+1][j+1] != -1 && hybrid_pos_target.substr(i+1,ches1[i+1][j+1]-i).find('1',0) == string::npos &&
           (Cseed[i+1][j+1]==MAX_DOUBLE ? MAX_DOUBLE : 
           (hybrid_energy=Cseed[i+1][j+1]+AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])))
			< min_hybrid_energy)
        {
         min_hybrid_energy = hybrid_energy;
         hybridization_found = true;
         hybrid_start_target = i+1;
         hybrid_start_miRNA = j+1;
         pos_target_1_dangle = i+1;
         pos_miRNA_1_dangle =j+1;
		}
	}
  }


  // mark positions of hybridization with '1' as hybridized in hybrid_pos_target
  // perform traceback
  if (hybridization_found)
   {
      traceback_print_labels_blanks(hybrid_start_target,hybrid_start_miRNA);
      pos_target_1=hybrid_start_target;
      pos_miRNA_1 =hybrid_start_miRNA;
      hybrid_pos_target.replace(pos_target_1_dangle, ches1[pos_target_1][pos_miRNA_1]-pos_target_1_dangle+1,
                                ches1[pos_target_1][pos_miRNA_1]-pos_target_1_dangle+1, '1');
      traceback_left(hybrid_start_target,hybrid_start_miRNA);
   }

  return pair<bool, double>(hybridization_found,min_hybrid_energy);
}

/****************************************************************************
* traceback in C starting in (i,j) (sequence characters for i and j are     *
* already written to hybrid_sequence_target and hybrid_sequence_miRNA in    *
* traceback_seed(int, int, int))                                            *
****************************************************************************/

void traceback_right(int i, int j)
{
 // cout << "traceback right (" << i << "," << j << ")" << endl;

  if ( i==(int)strlen(sequence_target) || j==(int)strlen(sequence_miRNA) || C[i][j]==MAX_DOUBLE )
    {
      cerr << "traceback_right error" << endl;
      cerr << "i=" << i << ",j=" << j << ",C=" << C[i][j] << endl;
      exit(1);
    }
  // last base-pair (i,j)
  else if ( i<(int)strlen(sequence_target)-1 && j<(int)strlen(sequence_miRNA)-1 && fabs(C[i][j] - 
             (ed5(int_sequence_miRNA[j],int_sequence_target[i],int_sequence_miRNA[j+1])+
              ed3(int_sequence_miRNA[j],int_sequence_target[i],int_sequence_target[i+1])+
              ED_miRNA[j][j+1]+ED_target[i][i+1]+AU_penalty(sequence_miRNA[j],sequence_target[i])+
              +duplex_init) ) <DOUBLE_DIFF ) {
         traceback_print_hybrid_end(i,j);
         pos_target_4=i;
         pos_miRNA_4 =j;
         pos_target_4_dangle =i+1;
         pos_miRNA_4_dangle =j+1;
       }
  else if ( i<(int)strlen(sequence_target)-1 && fabs(C[i][j] - 
             (ed3(int_sequence_miRNA[j],int_sequence_target[i],int_sequence_target[i+1])+
              ED_miRNA[j][j]+ED_target[i][i+1]+AU_penalty(sequence_miRNA[j],sequence_target[i])+
              +duplex_init) ) <DOUBLE_DIFF ) {
         traceback_print_hybrid_end(i,j);
         pos_target_4=i;
         pos_miRNA_4 =j;
         pos_target_4_dangle =i+1;
         pos_miRNA_4_dangle =j;
       }
  else if ( j<(int)strlen(sequence_miRNA)-1 && fabs(C[i][j] - 
             (ed5(int_sequence_miRNA[j],int_sequence_target[i],int_sequence_miRNA[j+1])+
              ED_miRNA[j][j+1]+ED_target[i][i]+AU_penalty(sequence_miRNA[j],sequence_target[i])+
              +duplex_init) ) <DOUBLE_DIFF ) {
         traceback_print_hybrid_end(i,j);
         pos_target_4=i;
         pos_miRNA_4 =j;
         pos_target_4_dangle =i;
         pos_miRNA_4_dangle =j+1;
       }
  else if ( fabs(C[i][j] - (AU_penalty(sequence_miRNA[j],sequence_target[i])+ED_target[i][i]+ED_miRNA[j][j]+
             +duplex_init) ) <DOUBLE_DIFF ) {
         traceback_print_hybrid_end(i,j);
         pos_target_4=i;
         pos_miRNA_4 =j;
         pos_target_4_dangle =i;
         pos_miRNA_4_dangle =j;
       }

  // stack
  else if ( che1[i+1][j+1] != -1 && che2[i+1][j+1] != -1 && fabs(C[i][j]-
		 (StackingEnergy(int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1], int_sequence_miRNA[j+1])+
		  C[i+1][j+1]+
		  ED_target[i][che1[i+1][j+1]]-ED_target[i+1][che1[i+1][j+1]]+
		  ED_miRNA [j][che2[i+1][j+1]]-ED_miRNA [j+1][che2[i+1][j+1]])) <DOUBLE_DIFF )
    {
      bulges_sequence_target+=" ";
      hybrid_sequence_target+=sequence_target[i];      
      hybrid_sequence_miRNA +=sequence_miRNA[j];      
      bulges_sequence_miRNA +=" ";
      traceback_right(i+1,j+1);
    }
  else
    {
      bool traceback_step=false;
      for(int k=1;!traceback_step && k<=Minimum(MAX_BULGE_SIZE,(int)strlen(sequence_target)-i-2);k++)      
	if ( che1[i+1+k][j+1] != -1 && che2[i+1+k][j+1] != -1 && fabs(C[i][j]-
		  (BulgeEnergy(k,int_sequence_target[i],int_sequence_miRNA[j],int_sequence_target[i+1+k],int_sequence_miRNA[j+1])+
		   C[i+1+k][j+1]+
		   ED_target[i][che1[i+1+k][j+1]]-ED_target[i+1+k][che1[i+1+k][j+1]]+
		   ED_miRNA [j][che2[i+1+k][j+1]]-ED_miRNA [j+1  ][che2[i+1+k][j+1]])) <DOUBLE_DIFF )
	  {
	    bulges_sequence_target+=" ";
	    hybrid_sequence_target+=sequence_target[i];      
	    hybrid_sequence_miRNA +=sequence_miRNA[j];      
	    bulges_sequence_miRNA +=" ";

	    for(int k_=i+1;k_<=i+k;k_++)
	      {
		bulges_sequence_target+=sequence_target[k_];
		hybrid_sequence_target+=" ";      
		hybrid_sequence_miRNA +=" ";      
		bulges_sequence_miRNA +=" ";
	      }
	    traceback_right(i+1+k,j+1);
	    traceback_step=true;
	  }
      for(int l=1;!traceback_step && l<=Minimum(MAX_BULGE_SIZE,(int)strlen(sequence_miRNA)-j-2);l++)      
	if ( che1[i+1][j+1+l] != -1 && che2[i+1][j+1+l] != -1 && fabs(C[i][j]-
		  (BulgeEnergy(l,int_sequence_target[i],int_sequence_miRNA[j],int_sequence_target[i+1],int_sequence_miRNA[j+1+l])+
		   C[i+1][j+1+l]+
		   ED_target[i][che1[i+1][j+1+l]]-ED_target[i+1   ][che1[i+1][j+1+l]]+
		   ED_miRNA [j][che2[i+1][j+1+l]]-ED_miRNA [j+1+l ][che2[i+1][j+1+l]])) <DOUBLE_DIFF )
	  {
	    bulges_sequence_target+=" ";
	    hybrid_sequence_target+=sequence_target[i];      
	    hybrid_sequence_miRNA +=sequence_miRNA[j];      
	    bulges_sequence_miRNA +=" ";

	    for(int l_=j+1;l_<=j+l;l_++)
	      {
		bulges_sequence_target+=" ";
		hybrid_sequence_target+=" ";      
		hybrid_sequence_miRNA +=" ";      
		bulges_sequence_miRNA +=sequence_miRNA[l_];
	      }
	    traceback_right(i+1,j+1+l);
	    traceback_step=true;
	  }
      for(int k=1;!traceback_step && k<=Minimum(MAX_INTERIOR_ONE_SIDE_SIZE,(int)strlen(sequence_target)-i-2);k++)
	for(int l=1;!traceback_step && l<=Minimum(MAX_INTERIOR_ONE_SIDE_SIZE,(int)strlen(sequence_miRNA)-j-2);l++)
	  if ( che1[i+1+k][j+1+l] != -1 && che2[i+1+k][j+1+l] != -1 && fabs(C[i][j]-
		    (InteriorLoopEnergy(k,l, int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1+k], int_sequence_miRNA[j+1+l], int_sequence_target[i+1], int_sequence_miRNA[j+1], int_sequence_target[i+k], int_sequence_miRNA[j+l])+ 
		     C[i+1+k][j+1+l]+
		     ED_target[i][che1[i+1+k][j+1+l]]-ED_target[i+1+k ][che1[i+1+k][j+1+l]]+
		     ED_miRNA [j][che2[i+1+k][j+1+l]]-ED_miRNA [j+1+l ][che2[i+1+k][j+1+l]])) <DOUBLE_DIFF )
	    {
	      bulges_sequence_target+=" ";
	      hybrid_sequence_target+=sequence_target[i];      
	      hybrid_sequence_miRNA +=sequence_miRNA[j];      
	      bulges_sequence_miRNA +=" ";

	      for(int k_=i+1;k_<=i+k;k_++)
		{
		  bulges_sequence_target+=sequence_target[k_];
		  hybrid_sequence_target+=" ";
		}
	      for(int k_=1;k_<=(k<l?l-k:0);k_++)
		{
		  bulges_sequence_target+=" ";
		  hybrid_sequence_target+=" ";
		}

	      for(int l_=j+1;l_<=j+l;l_++)
		{
		  hybrid_sequence_miRNA+=" ";
		  bulges_sequence_miRNA+=sequence_miRNA[l_];
		}
	      for(int l_=1;l_<=(l<k?k-l:0);l_++)
		{
		  hybrid_sequence_miRNA+=" ";
		  bulges_sequence_miRNA+=" ";;
		}

	      traceback_right(i+1+k,j+1+l);
	      traceback_step=true;
	    }     

      if (!traceback_step)
	{
	  cerr << "traceback_right error" << endl;
	  cerr << "i=" << i << ",j=" << j << ",C=" << C[i][j] << endl;
	  exit(1);
	}
    }
}

/****************************************************************************
* traceback in Cseed (matrix) starting in (i,j)                                 *
****************************************************************************/

void traceback_left(int i, int j)
{
  //  cout << "traceback left (" << i << "," << j << ")" << endl;

  if ( i==(int)strlen(sequence_target) || j==(int)strlen(sequence_miRNA) || Cseed[i][j]==MAX_DOUBLE )
    {
      cerr << "traceback_left error" << endl;
      cerr << "i=" << i << ",j=" << j << ",Cseed=" << Cseed[i][j] << endl;
      exit(1);
    }
  else if ( ches1[i+1][j+1] != -1 && ches2[i+1][j+1] != -1 && fabs(Cseed[i][j]-
		 (StackingEnergy(int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1], int_sequence_miRNA[j+1])+
		  Cseed[i+1][j+1]+
		  ED_target[i][ches1[i+1][j+1]]-ED_target[i+1][ches1[i+1][j+1]]+
		  ED_miRNA [j][ches2[i+1][j+1]]-ED_miRNA [j+1][ches2[i+1][j+1]])) <DOUBLE_DIFF )
    {
      bulges_sequence_target+=" ";
      hybrid_sequence_target+=sequence_target[i];      
      hybrid_sequence_miRNA +=sequence_miRNA[j];      
      bulges_sequence_miRNA +=" ";
      traceback_left(i+1,j+1);
    }
  else
    {
      bool traceback_step=false;
      for(int k=1;!traceback_step && k<=Minimum(MAX_BULGE_SIZE,(int)strlen(sequence_target)-i-2);k++)      
	if ( ches1[i+1+k][j+1] != -1 && ches2[i+1+k][j+1] != -1 && fabs(Cseed[i][j]-
		  (BulgeEnergy(k,int_sequence_target[i],int_sequence_miRNA[j],int_sequence_target[i+1+k],int_sequence_miRNA[j+1])+
		   Cseed[i+1+k][j+1]+
		   ED_target[i][ches1[i+1+k][j+1]]-ED_target[i+1+k][ches1[i+1+k][j+1]]+
		   ED_miRNA [j][ches2[i+1+k][j+1]]-ED_miRNA [j+1  ][ches2[i+1+k][j+1]])) <DOUBLE_DIFF )
	  {
	    bulges_sequence_target+=" ";
	    hybrid_sequence_target+=sequence_target[i];      
	    hybrid_sequence_miRNA +=sequence_miRNA[j];      
	    bulges_sequence_miRNA +=" ";

	    for(int k_=i+1;k_<=i+k;k_++)
	      {
		bulges_sequence_target+=sequence_target[k_];
		hybrid_sequence_target+=" ";      
		hybrid_sequence_miRNA +=" ";      
		bulges_sequence_miRNA +=" ";
	      }
	    traceback_left(i+1+k,j+1);
	    traceback_step=true;
	  }
      for(int l=1;!traceback_step && l<=Minimum(MAX_BULGE_SIZE,(int)strlen(sequence_miRNA)-j-2);l++)      
	if ( ches1[i+1][j+1+l] != -1 && ches2[i+1][j+1+l] != -1 && fabs(Cseed[i][j]-
		  (BulgeEnergy(l,int_sequence_target[i],int_sequence_miRNA[j],int_sequence_target[i+1],int_sequence_miRNA[j+1+l])+
		   Cseed[i+1][j+1+l]+
		   ED_target[i][ches1[i+1][j+1+l]]-ED_target[i+1   ][ches1[i+1][j+1+l]]+
		   ED_miRNA [j][ches2[i+1][j+1+l]]-ED_miRNA [j+1+l ][ches2[i+1][j+1+l]])) <DOUBLE_DIFF )
	  {
	    bulges_sequence_target+=" ";
	    hybrid_sequence_target+=sequence_target[i];      
	    hybrid_sequence_miRNA +=sequence_miRNA[j];      
	    bulges_sequence_miRNA +=" ";

	    for(int l_=j+1;l_<=j+l;l_++)
	      {
		bulges_sequence_target+=" ";
		hybrid_sequence_target+=" ";      
		hybrid_sequence_miRNA +=" ";      
		bulges_sequence_miRNA +=sequence_miRNA[l_];
	      }
	    traceback_left(i+1,j+1+l);
	    traceback_step=true;
	  }
      for(int k=1;!traceback_step && k<=Minimum(MAX_INTERIOR_ONE_SIDE_SIZE,(int)strlen(sequence_target)-i-2);k++)
	for(int l=1;!traceback_step && l<=Minimum(MAX_INTERIOR_ONE_SIDE_SIZE,(int)strlen(sequence_miRNA)-j-2);l++)
	  if ( ches1[i+1+k][j+1+l] != -1 && ches2[i+1+k][j+1+l] != -1 && fabs(Cseed[i][j]-
		    (InteriorLoopEnergy(k,l, int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1+k], int_sequence_miRNA[j+1+l], int_sequence_target[i+1], int_sequence_miRNA[j+1], int_sequence_target[i+k], int_sequence_miRNA[j+l])+ 
		     Cseed[i+1+k][j+1+l]+
		     ED_target[i][ches1[i+1+k][j+1+l]]-ED_target[i+1+k ][ches1[i+1+k][j+1+l]]+
		     ED_miRNA [j][ches2[i+1+k][j+1+l]]-ED_miRNA [j+1+l ][ches2[i+1+k][j+1+l]])) <DOUBLE_DIFF )
	    {
	      bulges_sequence_target+=" ";
	      hybrid_sequence_target+=sequence_target[i];      
	      hybrid_sequence_miRNA +=sequence_miRNA[j];      
	      bulges_sequence_miRNA +=" ";

	      for(int k_=i+1;k_<=i+k;k_++)
		{
		  bulges_sequence_target+=sequence_target[k_];
		  hybrid_sequence_target+=" ";
		}
	      for(int k_=1;k_<=(k<l?l-k:0);k_++)
		{
		  bulges_sequence_target+=" ";
		  hybrid_sequence_target+=" ";
		}

	      for(int l_=j+1;l_<=j+l;l_++)
		{
		  hybrid_sequence_miRNA+=" ";
		  bulges_sequence_miRNA+=sequence_miRNA[l_];
		}
	      for(int l_=1;l_<=(l<k?k-l:0);l_++)
		{
		  hybrid_sequence_miRNA+=" ";
		  bulges_sequence_miRNA+=" ";
		}

	      traceback_left(i+1+k,j+1+l);
	      traceback_step=true;
	    }     
      if (!traceback_step)
	{
	  if (!PairTest(sequence_target[i],sequence_miRNA[j]))
	    {
	      cerr << "traceback_left error" << endl;
	      cerr << "i=" << i << ",j=" << j << ": cannot pair" << endl;
	      exit(1);
	    }
	  
	  /*traceback_seed wants the bases x_i and y_j to be done, but 
	    traceback_left usually does not do this. thus, we have to 
	    output them here.*/
	  bulges_sequence_target+=" ";
	  hybrid_sequence_target+=sequence_target[i];      
	  hybrid_sequence_miRNA +=sequence_miRNA[j];
	  bulges_sequence_miRNA +=" ";

	  pos_target_2=i;
	  pos_miRNA_2 =j;

	  traceback_seed(i,j);
	}
    }
}

/****************************************************************************
* traceback in U (matrix) starting in (i,j)                                 *
****************************************************************************/

void traceback_U(int i, int j)
{
  //  cout << "traceback U (" << i << "," << j << ")" << endl;

  if ( i>=(int)strlen(sequence_target)-1 || j>=(int)strlen(sequence_miRNA)-1 )
    {
      cerr << "i(=" << i << ") or j(=" << j << ") out of range" << endl;
      exit(1);
    }

  if (fabs(U[i][j]-U[i+1][j])<DOUBLE_DIFF)
    traceback_U(i+1,j);
  else if (fabs(U[i][j]-U[i][j+1])<DOUBLE_DIFF)
    traceback_U(i,j+1);
  else if (fabs(U[i][j]-
            (Cseed[i+1][j+1]+ed5(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_target[i])+ed3(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_miRNA[j])
            +AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])+ED_target[i][ches1[i+1][j+1]]+ED_miRNA[j][ches2[i+1][j+1]]
            -ED_target[i+1][ches1[i+1][j+1]]-ED_miRNA[j+1][ches2[i+1][j+1]]))<DOUBLE_DIFF)
         {
            traceback_print_labels_blanks(i+1,j+1);
            pos_target_1=i+1;
            pos_miRNA_1 =j+1;
            pos_target_1_dangle = i;
            pos_miRNA_1_dangle =j;
            traceback_left(i+1,j+1);
         }
   else if (fabs(U[i][j]-
             (Cseed[i+1][j+1]+ed5(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_target[i])+AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])
             +ED_target[i][ches1[i+1][j+1]]-ED_target[i+1][ches1[i+1][j+1]]))<DOUBLE_DIFF)
         {
            traceback_print_labels_blanks(i+1,j+1);
            pos_target_1=i+1;
            pos_miRNA_1 =j+1;
            pos_target_1_dangle = i;
            pos_miRNA_1_dangle =j+1;
            traceback_left(i+1,j+1);
         }
   else if (fabs(U[i][j]-
             (Cseed[i+1][j+1]+ed3(int_sequence_target[i+1],int_sequence_miRNA[j+1],int_sequence_miRNA[j])+AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])
             +ED_miRNA[j][ches2[i+1][j+1]]-ED_miRNA[j+1][ches2[i+1][j+1]]))<DOUBLE_DIFF)
         {
            traceback_print_labels_blanks(i+1,j+1);
            pos_target_1=i+1;
            pos_miRNA_1 =j+1;
            pos_target_1_dangle = i+1;
            pos_miRNA_1_dangle =j;
            traceback_left(i+1,j+1);
         }
   else if (fabs(U[i][j]-(Cseed[i+1][j+1]+AU_penalty(sequence_target[i+1],sequence_miRNA[j+1])))<DOUBLE_DIFF)
         {
            traceback_print_labels_blanks(i+1,j+1);
            pos_target_1=i+1;
            pos_miRNA_1 =j+1;
            pos_target_1_dangle = i+1;
            pos_miRNA_1_dangle =j+1;
            traceback_left(i+1,j+1);
      }

  else
    {
      cerr << "traceback_U error" << endl;
      exit(1);
    }
}


/****************************************************************************
* traceback in hybrid_matrix starting in (i,j) (sequence characters for i   *
* and j are already written to hybrid_sequence_target and                   *
* hybrid_sequence_miRNA in traceback_seed(int, int))                        *
****************************************************************************/

void traceback_hybrid(int i, int j, int num_bp, int unpaired_target, int unpaired_miRNA)
{
  //  cout << "traceback hybrid (" << i << "," << j << ")" << endl;

  double stacking, bulge_target, bulge_miRNA, interior_loop;
  double energy = hybrid_matrix[i][j][num_bp][unpaired_target][unpaired_miRNA];
   
   if (num_bp == 2)
   {   
      /*stacking*/
      if ((unpaired_target == 0) && (unpaired_miRNA == 0) && (fabs(energy-StackingEnergy(int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1], int_sequence_miRNA[j+1])) < DOUBLE_DIFF))
      {
         bulges_sequence_target += " ";
         hybrid_sequence_target += sequence_target[i+1]; 
         hybrid_sequence_miRNA  += sequence_miRNA[j+1];
         bulges_sequence_miRNA  += " ";
      }            
               
      /*target_bulge*/         
      else if ((unpaired_target > 0) && (unpaired_target <= MAX_BULGE_SIZE) && (unpaired_miRNA == 0) && (fabs(energy-BulgeEnergy(unpaired_target, int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1+unpaired_target], int_sequence_miRNA[j+1])) < DOUBLE_DIFF))
      {         
         /*  -   x_{i+1} x_{i+2} ..... x_{i+unpaired_target} - 
            x_i    -        -    .....       -               x_{i+unpaired_target+1}
            y_j    -        -    .....       -               y_{j+1}
             -     -        -    .....       -               - */
             
         for (int h=1; h<=unpaired_target; h++)
         {
            bulges_sequence_target += sequence_target[i+h];
            hybrid_sequence_target += " "; 
            hybrid_sequence_miRNA  += " ";
            bulges_sequence_miRNA  += " ";
         }
         bulges_sequence_target += " ";
         hybrid_sequence_target += sequence_target[i+unpaired_target+1]; 
         hybrid_sequence_miRNA  += sequence_miRNA[j+1];
         bulges_sequence_miRNA  += " ";
      }         
           
      /*miRNA_bulge*/
      else if ((unpaired_target == 0) && (unpaired_miRNA <= MAX_BULGE_SIZE) && (unpaired_miRNA > 0) && (fabs(energy-BulgeEnergy(unpaired_miRNA, int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1], int_sequence_miRNA[j+1+unpaired_miRNA])) < DOUBLE_DIFF))
      {
         /*  -     -        -   .....       -                - 
            x_i    -        -   .....       -               x_{i+1}
            y_j    -        -   .....       -               y_{j+unpaired_miRNA+1}
             -  y_{j+1} y_{j+2} ..... y_{j+unpaired_miRNA}   - */
             
         for (int h=1; h<=unpaired_miRNA; h++)
         {
            bulges_sequence_target += " ";
            hybrid_sequence_target += " "; 
            hybrid_sequence_miRNA  += " ";
            bulges_sequence_miRNA  += sequence_miRNA[j+h];
         }
         bulges_sequence_target += " ";
         hybrid_sequence_target += sequence_target[i+1]; 
         hybrid_sequence_miRNA  += sequence_miRNA[j+unpaired_miRNA+1];
         bulges_sequence_miRNA  += " ";
      }
                  
      /* interior loop */
      else if ((unpaired_target > 0) && (unpaired_target <= MAX_INTERIOR_ONE_SIDE_SIZE) && (unpaired_miRNA > 0) && (unpaired_miRNA <= MAX_INTERIOR_ONE_SIDE_SIZE) && (fabs(energy-InteriorLoopEnergy(unpaired_target, unpaired_miRNA, int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1+unpaired_target], int_sequence_miRNA[j+1+unpaired_miRNA], int_sequence_target[i+1], int_sequence_miRNA[j+1], int_sequence_target[i+unpaired_target], int_sequence_miRNA[j+unpaired_miRNA])) < DOUBLE_DIFF))
      {
      /*     -   x_{i+1} x_{i+2} .......................... x_{i+unpaired_target} - 
            x_i    -        -    ..........................       -               x_{i+unpaired_target+1}
            y_j    -        -    ..........................       -               y_{j+unpaired_miRNA+1}
             -   y_{j+1} y_{j+2} ... y_{j+unpaired_miRNA}         -               - */
             
         int bulge_max = Maximum(unpaired_target,unpaired_miRNA);
         for (int h=1; h<=bulge_max; h++)
         {
            hybrid_sequence_target += " "; 
            hybrid_sequence_miRNA  += " ";
            
            /* two cases since the interior loop might be asymmetric*/
            if (unpaired_target < h)
               bulges_sequence_target += " ";
            else
               bulges_sequence_target += sequence_target[i+h];
            
            if (unpaired_miRNA < h)
               bulges_sequence_miRNA += " ";
            else
               bulges_sequence_miRNA  += sequence_miRNA[j+h];
         }    
             
         bulges_sequence_target += " ";
         hybrid_sequence_target += sequence_target[i+unpaired_target+1]; 
         hybrid_sequence_miRNA  += sequence_miRNA[j+unpaired_miRNA+1];
         bulges_sequence_miRNA  += " ";  
         
      }
      
      else
      {
         cerr << "Error in the base case traceback of hybrid!" << endl;
      }
            
   }
   else if (num_bp > 2)
   {
      /*stacking*/
      stacking = Sum_or_MAX_DOUBLE(StackingEnergy(int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1], int_sequence_miRNA[j+1]), hybrid_matrix[i+1][j+1][num_bp-1][unpaired_target][unpaired_miRNA]);
      
      if (fabs(energy - stacking) < DOUBLE_DIFF)
      {
         bulges_sequence_target += " ";
         hybrid_sequence_target += sequence_target[i+1];
         hybrid_sequence_miRNA  += sequence_miRNA[j+1];
         bulges_sequence_miRNA  += " ";
         traceback_hybrid(i+1, j+1, num_bp-1, unpaired_target, unpaired_miRNA);
         return;
      }
         
       
      /* bulge in the target mRNA */         
      for (int p=1; p<=Minimum(unpaired_target,MAX_BULGE_SIZE); p++)
      {
         bulge_target = Sum_or_MAX_DOUBLE(BulgeEnergy(p, int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1+p], int_sequence_miRNA[j+1]), hybrid_matrix[i+1+p][j+1][num_bp-1][unpaired_target-p][unpaired_miRNA]);
            
         if (fabs(energy - bulge_target) < DOUBLE_DIFF)
         {
            /*  -   x_{i+1} x_{i+2} ..... x_{i+unpaired_target} - 
            x_i    -        -    .....       -               x_{i+unpaired_target+1}
            y_j    -        -    .....       -               y_{j+1}
             -     -        -    .....       -               - */
             
            for (int h=1; h<=p; h++)
            {
               bulges_sequence_target += sequence_target[i+h];
               hybrid_sequence_target += " "; 
               hybrid_sequence_miRNA  += " ";
               bulges_sequence_miRNA  += " ";
            }
            bulges_sequence_target += " ";
            hybrid_sequence_target += sequence_target[i+p+1]; 
            hybrid_sequence_miRNA  += sequence_miRNA[j+1];
            bulges_sequence_miRNA  += " ";
            
            traceback_hybrid(i+p+1,j+1,num_bp-1,unpaired_target-p, unpaired_miRNA);
            return;
         }

      }
         
      /* bulge in the miRNA */
      for (int q=1; q<=Minimum(unpaired_miRNA,MAX_BULGE_SIZE); q++)
      {
         bulge_miRNA = Sum_or_MAX_DOUBLE(BulgeEnergy(q, int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1], int_sequence_miRNA[j+1+q]), hybrid_matrix[i+1][j+1+q][num_bp-1][unpaired_target][unpaired_miRNA-q]);
            
         if (fabs(energy - bulge_miRNA) < DOUBLE_DIFF)
         {
            /*  -     -        -   .....       -                - 
            x_i    -        -   .....       -               x_{i+1}
            y_j    -        -   .....       -               y_{j+unpaired_miRNA+1}
             -  y_{j+1} y_{j+2} ..... y_{j+unpaired_miRNA}   - */
             
            for (int h=1; h<=q; h++)
            {
               bulges_sequence_target += " ";
               hybrid_sequence_target += " "; 
               hybrid_sequence_miRNA  += " ";
               bulges_sequence_miRNA  += sequence_miRNA[j+h];
            }
            bulges_sequence_target += " ";
            hybrid_sequence_target += sequence_target[i+1]; 
            hybrid_sequence_miRNA  += sequence_miRNA[j+q+1];
            bulges_sequence_miRNA  += " ";
            
            traceback_hybrid(i+1,j+q+1,num_bp-1,unpaired_target, unpaired_miRNA-q);
            return;
         }
      }
         
      /* interior loop */
      for (int p=1; p<=Minimum(unpaired_target,MAX_INTERIOR_ONE_SIDE_SIZE); p++)
      {
         for (int q=1; q<=Minimum(unpaired_miRNA,MAX_INTERIOR_ONE_SIDE_SIZE); q++)
         {
            interior_loop = Sum_or_MAX_DOUBLE(InteriorLoopEnergy(p, q, int_sequence_target[i], int_sequence_miRNA[j], int_sequence_target[i+1+p], int_sequence_miRNA[j+1+q], int_sequence_target[i+1], int_sequence_miRNA[j+1], int_sequence_target[i+p], int_sequence_miRNA[j+q]), hybrid_matrix[i+1+p][j+1+q][num_bp-1][unpaired_target-p][unpaired_miRNA-q]);
               
            if (fabs(energy - interior_loop) < DOUBLE_DIFF)
            {
               /*     -   x_{i+1} x_{i+2} .......................... x_{i+unpaired_target} - 
               x_i    -        -    ..........................       -               x_{i+unpaired_target+1}
               y_j    -        -    ..........................       -               y_{j+unpaired_miRNA+1}
                -   y_{j+1} y_{j+2} ... y_{j+unpaired_miRNA}         -               - */
             
               for (int h=1; h<=Maximum(p,q); h++)
               {
                  hybrid_sequence_target += " "; 
                  hybrid_sequence_miRNA  += " ";
            
                  /* two cases since the interior loop might be asymmetric*/
                  if (p < h)
                     bulges_sequence_target += " ";
                  else
                     bulges_sequence_target += sequence_target[i+h];
            
                  if (q < h)
                     bulges_sequence_miRNA += " ";
                  else
                     bulges_sequence_miRNA += sequence_miRNA[j+h];
               }    
             
               bulges_sequence_target += " ";
               hybrid_sequence_target += sequence_target[i+p+1]; 
               hybrid_sequence_miRNA  += sequence_miRNA[j+q+1];
               bulges_sequence_miRNA  += " ";  
               
               traceback_hybrid(i+p+1,j+q+1, num_bp-1, unpaired_target-p, unpaired_miRNA-q);
               return;
            }
         }
      }
         
   }
   else
      cerr << "Error in hybrid traceback!" << endl;
}



/****************************************************************************
* traceback in seed starting in (i,j) (sequence characters for i and j are  *
* already written to hybrid_sequence_target and hybrid_sequence_miRNA in    *
* traceback(int))                                                           *
****************************************************************************/

void traceback_seed(int i, int j)
{
   if (seed[i][j] == MAX_DOUBLE)
   {
      cerr << "No seed found!\n";
   }
   else
   {
      //traceback of hybrid
      traceback_hybrid( i, j, num_bp_seed, seed_unpaired_target[i][j], seed_unpaired_miRNA[i][j]);
      //traceback of C

      /* the following has to be done, since the traceback of hybrid has already analyzed bases
         x_{i+num_bp_seed-1+seed_unpaired_target[i][j]} and y_{j+num_bp_seed-1+seed_unpaired_miRNA[i][j]} 
         but traceback_right also does it. therefore, we have to remove the last bases here.*/
      bulges_sequence_target = bulges_sequence_target.substr(0,bulges_sequence_target.length()-1);
      hybrid_sequence_target = hybrid_sequence_target.substr(0,hybrid_sequence_target.length()-1);
      hybrid_sequence_miRNA  = hybrid_sequence_miRNA .substr(0,hybrid_sequence_miRNA .length()-1); 
      bulges_sequence_miRNA  = bulges_sequence_miRNA .substr(0,bulges_sequence_miRNA .length()-1);  

      pos_target_3=i+num_bp_seed-1+seed_unpaired_target[i][j];
      pos_miRNA_3 =j+num_bp_seed-1+seed_unpaired_miRNA[i][j] ;
      
      traceback_right( i+num_bp_seed-1+seed_unpaired_target[i][j], j+num_bp_seed-1+seed_unpaired_miRNA[i][j]);
   }
   
   return;
}

/****************************************************************************
* traceback in Cseed (matrix) starting in (0,0)                             *
****************************************************************************/

bool traceback(double energy)
{
  // recompute matrices for best che value before traceback
  if (!heuristic) {
    init_global_matrices();

    for(int i=0;i<=best_che1;i++)
      for(int j=0;j<=best_che2;j++) {
        che1[i][j] = best_che1;
        che2[i][j] = best_che2;
      }

    // compute matrices for all indices <= best_che1,best_che2
    for(int i=best_che1;i>=0;i--)
      for(int j=best_che2;j>=0;j--)
        C[i][j]=Cij(i,j);

    for (int i=best_che1;i>=0;i--)
      for(int j=best_che2;j>=0;j--)
        for (int l=2; l<=num_bp_seed; l++)
          for (int b_target=0; b_target<=max_unpaired_target; b_target++)
            for (int b_miRNA=0; b_miRNA<=max_unpaired_miRNA; b_miRNA++)
              hybrid_matrix[i][j][l][b_target][b_miRNA] = calc_hybrid_value(i,j,l,b_target,b_miRNA);

    for(int i=best_che1;i>=0;i--)
      for(int j=best_che2;j>=0;j--) {
        seed [i][j] = compute_seed_ij(i,j);
        Cseed[i][j] = Cij_seed(i,j);
        U    [i][j] = U_ij(i,j);
      }
  }


  if (fabs(energy-0.0)<DOUBLE_DIFF)
    {
      return false;
    }

  else {
   bool traceback_step=false;

   if (fabs(energy-(Cseed[0][0]+AU_penalty(sequence_target[0],sequence_miRNA[0])))<DOUBLE_DIFF) {
      traceback_print_labels_blanks(0,0);
      pos_target_1=0;
      pos_miRNA_1 =0;
      pos_target_1_dangle = 0;
      pos_miRNA_1_dangle =0;
      traceback_left(0,0);

      traceback_step = true;
   }

   for (int i=1; !traceback_step && i<(int)strlen(sequence_target);i++)
      {
        if (fabs(energy-(Cseed[i][0]+ed5(int_sequence_target[i],int_sequence_miRNA[0],int_sequence_target[i-1])+AU_penalty(sequence_target[i],sequence_miRNA[0])
                         +ED_target[i-1][ches1[i][0]]-ED_target[i][ches1[i][0]]))<DOUBLE_DIFF)
          {
            traceback_print_labels_blanks(i,0);

            pos_target_1=i;
            pos_miRNA_1 =0;
            pos_target_1_dangle = i-1;
            pos_miRNA_1_dangle =0;
            traceback_left(i,0);

            traceback_step = true;
          }
        else if (fabs(energy-(Cseed[i][0]+AU_penalty(sequence_target[i],sequence_miRNA[0])))<DOUBLE_DIFF)
          {
            traceback_print_labels_blanks(i,0);

            pos_target_1=i;
            pos_miRNA_1 =0;
            pos_target_1_dangle = i;
            pos_miRNA_1_dangle =0;
            traceback_left(i,0);

            traceback_step = true;
          }
      }

     for (int j=1; !traceback_step && j<(int)strlen(sequence_miRNA);j++)
      {
        if (fabs(energy-(Cseed[0][j]+ed3(int_sequence_target[0],int_sequence_miRNA[j],int_sequence_miRNA[j-1])+AU_penalty(sequence_target[0],sequence_miRNA[j])
                         +ED_miRNA[j-1][ches2[0][j]]-ED_miRNA[j][ches2[0][j]]))<DOUBLE_DIFF)
          {
            traceback_print_labels_blanks(0,j);

            pos_target_1=0;
            pos_miRNA_1 =j;
            pos_target_1_dangle = 0;
            pos_miRNA_1_dangle =j-1;
            traceback_left(0,j);

            traceback_step = true;
          }
        else if (fabs(energy-(Cseed[0][j]+AU_penalty(sequence_target[0],sequence_miRNA[j])))<DOUBLE_DIFF)
          {
            traceback_print_labels_blanks(0,j);

            pos_target_1=0;
            pos_miRNA_1 =j;
            pos_target_1_dangle = 0;
            pos_miRNA_1_dangle =j;
            traceback_left(0,j);

            traceback_step = true;
          }

      }

      if (!traceback_step)
	{
          traceback_U(0,0);
        }
   }

   // init hybrid_pos_target to allow search for suboptimal hybridizations
   hybrid_pos_target.replace(pos_target_1_dangle, ches1[pos_target_1][pos_miRNA_1]-pos_target_1_dangle+1,
                             ches1[pos_target_1][pos_miRNA_1]-pos_target_1_dangle+1, '1');

  return true;
}

/****************************************************************************
* add the direction label and insert blanks at beginning of the             *
* hybridization sequence and the bulge sequence                             *
* i,j are hybridization start in target, miRNA
****************************************************************************/

void traceback_print_labels_blanks(int i, int j)
{
  // insert spaces in target string
  if (i<j)
    {
      for(int k=0;k<j-i;k++)
       {
         bulges_sequence_target+=" ";
         hybrid_sequence_target+=" ";
       }
       bulges_sequence_target+="5'-";
       hybrid_sequence_target+="   ";
       for(int k=0;k<i;k++)
        {
          bulges_sequence_target+=sequence_target[k];
          hybrid_sequence_target+=" ";
        }
        hybrid_sequence_miRNA +="   ";
        bulges_sequence_miRNA +="3'-";
        for(int l=0;l<j;l++)
         {
           hybrid_sequence_miRNA+=" ";
           bulges_sequence_miRNA+=sequence_miRNA[l];
         }
    }
   // insert spaces in miRNA string
   else
      {
        for(int l=0;l<i-j;l++)
         {
           hybrid_sequence_miRNA+=" ";      
           bulges_sequence_miRNA+=" ";
         }
        hybrid_sequence_miRNA+="   ";
        bulges_sequence_miRNA+="3'-";
        for(int l=0;l<j;l++)
         {
           hybrid_sequence_miRNA+=" ";      
           bulges_sequence_miRNA+=sequence_miRNA[l];
         }
        bulges_sequence_target +="5'-";;
        hybrid_sequence_target +="   ";      
        for(int k=0;k<i;k++)
         {
           bulges_sequence_target+=sequence_target[k];
           hybrid_sequence_target+=" ";      
         }
      }
}

/****************************************************************************
* prints last right base pair of hybridization, prints remaining non-       *
* hybridized sequence and adds direction labels                             *
****************************************************************************/

void traceback_print_hybrid_end(int i, int j)
{
      bulges_sequence_target+=" ";
      hybrid_sequence_target+=sequence_target[i];
      hybrid_sequence_miRNA +=sequence_miRNA[j];
      bulges_sequence_miRNA +=" ";
 
      for(int k=i+1;k<(int)strlen(sequence_target);k++)
	bulges_sequence_target+=sequence_target[k];
      for(int l=j+1;l<(int)strlen(sequence_miRNA);l++)
	bulges_sequence_miRNA +=sequence_miRNA[l];
      bulges_sequence_target+="-3'";
      bulges_sequence_miRNA +="-5'";
}
