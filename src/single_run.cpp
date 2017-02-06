#include "single_run.h"


/*************************************************************************************
* reads a fasta file of sequences in and stores them in a global vector of strings,  *
* the sequences are written to every second line of the vector, the other lines are  *
* filled with the sequences' names                                                   *
*************************************************************************************/


bool readin_fasta(char* filename, vector<string>& vec)
{
   bool okay;
   ifstream seqin(filename);
   
   if (seqin.fail())
   {
      cerr << "The file " << filename << " does not exist!\n";
      exit(1);
   }
   
   string tmp, tmp_seq = "";
   while (!seqin.eof())
   {
      getline(seqin,tmp);
      if (tmp[0] == '>')
      {
         if (tmp_seq == "")
	     {
	        if ((int)vec.size() == 0)
	           vec.push_back(tmp);
	        else // sequence missing
	        {
              cout << "# WARNING: Sequence " << vec[(int)vec.size()-1] << " is missing and thus skipped.\n";
       	       vec.pop_back();
                vec.push_back(tmp);
	        }
	     }
	     else 
	     {   
	        okay = check_sequence_and_set2upper(tmp_seq);
	        if (okay)
	           vec.push_back(tmp_seq);
	        else
            {
              cout << "# WARNING: Sequence " << vec[(int)vec.size()-1] << " in file " << filename << " is skipped due to invalid character(s).\n";
	           vec.pop_back();			   
	        }
	        vec.push_back(tmp);
	        tmp_seq = "";
	     }
      }
      else
	     tmp_seq += tmp;
   }
   
   seqin.close();
   
   //if file empty
   if ((int)vec.size() == 0)
      return false;

   
   //final sequence
   if (tmp_seq == "")
   {
      cout << "# WARNING: Sequence " << vec[(int)vec.size()-1] << " is missing and thus skipped.\n";
      vec.pop_back();
   }
   else
   {
      okay = check_sequence_and_set2upper(tmp_seq);
      if (okay)
         vec.push_back(tmp_seq);
      else
      {
     cout << "# WARNING: Sequence " << vec[(int)vec.size()-1] << " in file " << filename << " is skipped due to invalid character(s).\n";
	 vec.pop_back();			   
      }
   }
   
   //print
   // for (int i=0; i<(int)vec.size(); i++)
   //   cout << vec[i] << endl;
   
	  
   // if size == 0 => no sequence in the file or all sequences have errors
   if ((int)vec.size() == 0)
      return false;
   else
      return true;
}



/***************************************************************************
* test, whether the values of the input options relevant for the           *
  hybridization are valid for the current sequences                        *
***************************************************************************/

bool test_options_hybridization()
{
   
   if (num_bp_seed < 2)
   {
      cerr << "The length of the seed has to be at least 2!\n\n";
      return false;
   }
   
   if (num_bp_seed > Minimum((int)strlen(sequence_target),(int)strlen(sequence_miRNA)))
   {
      cerr << "The length of the seed has to be less than or equal the length of the shortest input sequence!\n\n";
      return false;
   }
      
   if ((max_unpaired_target < 0) || (max_unpaired_miRNA < 0) || (max_unpaired < 0))
   {
      if (max_unpaired < 0)
         cerr << "The maximal sum of unpaired bases in the seed region has to be greater or equal to 0!\n\n";
      else
         cerr << "The maximal number of unpaired bases in the seed region of both sequences has to be greater or equal to 0!\n\n";
      return false;
   }

   if (max_unpaired_target > (int)strlen(sequence_target))
   {
      cerr << "The number of unpaired bases in the seed region of sequence 1 has to be less or equal to the length of sequence 1!\n\n";
      return false;
   }
   
   if (max_unpaired_miRNA > (int)strlen(sequence_miRNA))
   {
      cerr << "The number of unpaired bases in the seed region of sequence 2 has to be less or equal to the length of sequence 2!\n\n";
      return false;
   }
   
   if (max_unpaired > ((int)strlen(sequence_target)+(int)strlen(sequence_miRNA)))
   {
      cerr << "The number of unpaired bases in the seed region of both sequences has to be less or equal to the sum of the lengths of the sequences!\n\n";
      return false;
   }

   if (seed_region_ncRNA && num_bp_seed > seed_region_end_ncRNA_5-seed_region_start_ncRNA_5+1)
   {
      cerr << "The length of the seed has to be at most the length of the seed region in the binding RNA!\n\n";
      return false;
   }

   if  (seed_region_ncRNA && seed_region_start_ncRNA_5 > seed_region_end_ncRNA_5 )
   {
      cerr << "The start of the seed region in the binding RNA has to be less than or equal to the end of the seed region!\n\n";
      return false;
   }

   if  (seed_region_ncRNA && seed_region_start_ncRNA_5 < 1 )
   {
      cerr << "The start of the seed region in the binding RNA has to be at least 1!\n\n";
      return false;
   }

   if ( seed_region_ncRNA && seed_region_end_ncRNA_5 > (int)strlen(sequence_miRNA))
   {
      cerr << "The end of the seed region in the binding RNA has to be less than or equal to the length of its sequence!\n\n";
      return false;
   }

   if (max_unpaired > max_unpaired_target + max_unpaired_miRNA)
   {
      if (max_unpaired_target == 0 && max_unpaired_miRNA == 0)
      {
         max_unpaired_target = max_unpaired_miRNA = max_unpaired;
      }
      else if (max_unpaired_miRNA == 0)
      {
         max_unpaired_miRNA = max_unpaired;
      }
      else if (max_unpaired_target == 0)
      {
         max_unpaired_target = max_unpaired;
      }
   }

   if (max_unpaired == 0 && (max_unpaired_target != 0 || max_unpaired_miRNA != 0))
   {
      cerr << "The number of unpaired bases in the seed region of both sequences should be at least equal to the number of unpaired bases in a single sequence!\n\n";
      return false;
   }

   return true;
  
}


/***************************************************************************
* test, whether the values of the input options relevant for the           *
  sliding window are valid for the current sequences                       *
***************************************************************************/

bool test_options_window()
{
   if ((window == true) && (window_size < 3))   // RNAplfold crashes if W is set smaller than 3
   {
      cerr << "The window size has to be higher or equal to 3!\n";
	  return false;
   }
   
   if (window_size > (int)strlen(sequence_target))
   {
      cerr << "The window size has to be at most the size of the target sequence!\n";
	  return false;
   }
   
   if ((max_length == true) && (max_unpaired_length < 1))
   {
      cerr << "The max. length of unpaired bases has to be higher or equal to 1!\n";
	  return false;
   }
   
   if (max_unpaired_length > window_size)
   {
      cerr << "The max. length of unpaired bases has to be lower or equal to the window size!\n";
	  return false;
   }
   
   if ((max_dist == true) && (max_pair_dist < 1))
   {
      cerr << "The max. distance of two paired bases has to be higher or equal to 1!\n";
	  return false;
   }   
   
   if (max_pair_dist > window_size)
   {
      cerr << "The max. distance of two paired bases has to be lower or equal to the window size!\n";
	  return false;
   }
   
/*   if (max_unpaired_length < max_unpaired)
   {
      cerr << "The max. length of unpaired bases has to be higher or equal to the number of unpaired base in the seed region of both sequences!\n";
          return false;
   }
*/   
   return true;
}


/**********************************************************************************
*         Init all data structures for target                                     *
**********************************************************************************/

void init_target()
{
   //translate sequences into integer
   int_sequence_target = char2int(sequence_target);

   /*************************************************************************
   *  if window and/or max_length are false, the values for window_size    *
   *  and max_unpaired_length are set to the maximal values (in the target)*
   *  if max_dist is false, the value of max_bp_dist is set to the window  *
   *  size                                                                 *
   ************************************************************************/

   if (window == false)
      window_size = (int)strlen(sequence_target);

   if (max_length == false)
       // max_unpaired_length has to be longer than length of miRNA sequence, 
       // since interacting region in target can contain long bulges
      max_unpaired_length = window_size;

   if (max_dist == false)
      max_pair_dist = window_size;

   /***************************************************************************
   * test, whether the values of the input options relevant for              *
   * hybridization are valid                                                 *
   ***************************************************************************/

   if (!test_options_window())
     exit(1);


   allocate_ED_target_matrix();  /* space allocation */
   init_ED_target_matrix()    ;  /* init */

   ED_target_computed  = false;

//has been moved to compute_seed_ij and is computed only when valid seed has been found
//    // calculate ED matrix for target
//    if (use_RNAplfold)
//       calc_ED_target_matrix_plfold(sequence_target);
//    else
//    {
//       if (window)
//          calc_ED_target_matrix_window();
//       else
//         calc_ED_target_matrix(); 
//    }
}


/**********************************************************************************
*         Init all data structures for miRNA                                      *
**********************************************************************************/

void init_miRNA()
{
   //translate sequences into integer
   int_sequence_miRNA  = char2int(sequence_miRNA);

   allocate_ED_miRNA_matrix();  /* space allocation */
   init_ED_miRNA_matrix()    ;  /* init */

   ED_miRNA_computed = false;

//has been moved to compute_seed_ij and is computed only when valid seed has been found
//    // calculate ED matrix for miRNA
//    calc_ED_miRNA_matrix();
}


/**********************************************************************************
*         Free all data structures for miRNA                                      *
**********************************************************************************/

void free_target()
{
   free_ED_target_matrix();
   free(int_sequence_target);
}


/**********************************************************************************
*         Free all data structures for miRNA                                      *
**********************************************************************************/

void free_miRNA()
{
   free_ED_miRNA_matrix();
   free(int_sequence_miRNA);
}


/**********************************************************************************
*         Clear sequences after hybridization                                     *
**********************************************************************************/

void clear_hydrization_seq()
{
   bulges_sequence_target.clear();
   hybrid_sequence_target.clear();
   hybrid_sequence_miRNA.clear();
   bulges_sequence_miRNA.clear();
}


/**********************************************************************************
*         Print the input and the results of the initializing step                *
**********************************************************************************/

void print_input()
{
  cout << "-------------------------"                                                               << endl;
  cout << "INPUT "                                                                                  << endl;
  cout << "-------------------------"                                                               << endl;
  cout << "number of base pairs in seed  : "                                 << num_bp_seed         << endl;
  cout << "max. number of unpaired bases in the seed region of seq. 1    : " << max_unpaired_target << endl;
  cout << "max. number of unpaired bases in the seed region of seq. 2    : " << max_unpaired_miRNA  << endl;
  cout << "max. number of unpaired bases in the seed region of both seq's: " << max_unpaired        << endl;
  if (seed_region_ncRNA)
	 //recalculated start and end positions since up to now we are using the reverse of the miRNA
     cout << "search region for seed in seq. 2                              : " 
		  << seed_region_start_ncRNA_5  << " -- "
		  << seed_region_end_ncRNA_5                                         << endl;
  cout << "RNAup used                                                    : " << (use_RNAup ? "true":"false") << endl;
  cout << "RNAplfold used                                                : " << (use_RNAplfold ? "true":"false") << endl;
  if (window || use_RNAplfold)
     cout << "sliding window size                                           : " << window_size      << endl;
  if (max_length || use_RNAplfold)
     cout << "max. length of unpaired region                                : " << max_unpaired_length << endl;
  if (max_dist || use_RNAplfold)
     cout << "max. distance of two paired bases                             : " << max_pair_dist    << endl;
  cout << "weight for ED values of target RNA in energy                  : " << ED_target_scaling << endl;
  cout << "weight for ED values of binding RNA in energy                 : " << ED_miRNA_scaling  << endl;
  if (fabs(pu_comb_threshold+1.0) > DOUBLE_DIFF)
     cout << "threshold for seed accessibility                              : " << pu_comb_threshold  << endl;
  cout << "temperature                                                   : " << temperature << " Celsius" << endl;
  cout << "max. number of subopt. results                                : " << max_subopt  << endl;
  cout << "Heuristic for hybridization end used                          : " << (heuristic ? "true":"false") << endl << endl;
  cout << "-------------------------"                                                               << endl;
  cout << "OUTPUT "                                                                                 << endl;
  cout << "-------------------------"                                                               << endl;

  return;
}


void print_matrices()
{
  cout << "ED_target" << endl;
  for(int i=0;i<(int)strlen(sequence_target);i++)
    {
      for(int j=0;j<(int)strlen(sequence_target);j++)
	cout << ED_target[i][j] << "\t";
      cout << endl;
    }
  cout << endl;

  cout << "ED_ncRNA" << endl;
  for(int i=0;i<(int)strlen(sequence_miRNA);i++)
    {
      for(int j=0;j<(int)strlen(sequence_miRNA);j++)
	cout << ED_miRNA[i][j] << "\t";
      cout << endl;
    }
  cout << endl;

  cout << "C" << endl;
  for(int i=0;i<(int)strlen(sequence_target);i++)
    {
      for(int j=0;j<(int)strlen(sequence_miRNA);j++)
	cout << C[i][j] << "\t";
      cout << endl;
    }
  cout << endl;

  cout << "C_seed" << endl;
  for(int i=0;i<(int)strlen(sequence_target);i++)
    {
      for(int j=0;j<(int)strlen(sequence_miRNA);j++)
	cout << Cseed[i][j] << "\t";
      cout << endl;
    }
  cout << endl;
  
  cout << "U" << endl;
  for(int i=0;i<(int)strlen(sequence_target);i++)
    {
      for(int j=0;j<(int)strlen(sequence_miRNA);j++)
	cout << U[i][j] << "\t";
      cout << endl;
    }
  cout << endl;

  cout << "seed" << endl;
  for(int i=0;i<(int)strlen(sequence_target);i++)
    {
      for(int j=0;j<(int)strlen(sequence_miRNA);j++)
	cout << seed[i][j] << "\t";
      cout << endl;
    }
  cout << endl;
  
  cout << "seed_unpaired_target" << endl;
  for (int i=0; i<(int)strlen(sequence_target); i++)
  {
     for (int j=0; j<(int)strlen(sequence_miRNA); j++)
        cout << seed_unpaired_target[i][j] << "\t";
     cout << endl;
  }
  cout << endl;
  
  cout << "seed_unpaired_ncRNA" << endl;
  for (int i=0; i<(int)strlen(sequence_target); i++)
  {
     for (int j=0; j<(int)strlen(sequence_miRNA); j++)
        cout << seed_unpaired_miRNA[i][j] << "\t";
     cout << endl;
  }
  cout << endl;

  cout << "che1" << endl;
  for(int i=0;i<(int)strlen(sequence_target);i++)
    {
      for(int j=0;j<(int)strlen(sequence_miRNA);j++)
	cout << che1[i][j] << "\t";
      cout << endl;
    }
  cout << endl;

  cout << "che2" << endl;
  for(int i=0;i<(int)strlen(sequence_target);i++)
    {
      for(int j=0;j<(int)strlen(sequence_miRNA);j++)
	cout << che2[i][j] << "\t";
      cout << endl;
    }
  cout << endl;

  cout << "ches1" << endl;
  for(int i=0;i<(int)strlen(sequence_target);i++)
    {
      for(int j=0;j<(int)strlen(sequence_miRNA);j++)
	cout << ches1[i][j] << "\t";
      cout << endl;
    }
  cout << endl;

  cout << "ches2" << endl;
  for(int i=0;i<(int)strlen(sequence_target);i++)
    {
      for(int j=0;j<(int)strlen(sequence_miRNA);j++)
	cout << ches2[i][j] << "\t";
      cout << endl;
    }
  cout << endl;
}


/**********************************************************************************
*         Print hybridization                                                     *
**********************************************************************************/

void print_hybridization(double v)
{
  cout << bulges_sequence_target               << endl;
  cout << hybrid_sequence_target               << endl;
  cout << hybrid_sequence_miRNA                << endl; 
  cout << bulges_sequence_miRNA        << endl << endl; 

  if (detailed_output)
    {
      //cout << "positions(target): -- " << pos_target_1 << " -- " << pos_target_2 << " -- " << pos_target_3 << " -- " << pos_target_4 << endl;
      //cout << "positions(miRNA ): -- " << pos_miRNA_1  << " -- " << pos_miRNA_2  << " -- " << pos_miRNA_3  << " -- " << pos_miRNA_4  << endl;

     //pos+1 for using a numbering starting at 1
     cout << "positions(target)     : " << pos_target_1+1 << " -- " << pos_target_4+1 << endl;
     cout << "positions seed(target): " << pos_target_2+1 << " -- " << pos_target_3+1 << endl;
     //positions including dangling ends
     cout << "positions with dangle(target): " << pos_target_1_dangle+1 << " -- " << pos_target_4_dangle+1 << endl;
     //recalculated start and end positions since up to now we are using the reverse of the miRNA
     int real_pos_miRNA_1 = (int)strlen(sequence_miRNA) - pos_miRNA_1;
     int real_pos_miRNA_4 = (int)strlen(sequence_miRNA) - pos_miRNA_4;
     cout << "positions(ncRNA)      : " << real_pos_miRNA_4  << " -- " << real_pos_miRNA_1  << endl;
     cout << "positions seed(ncRNA) : " << (int)strlen(sequence_miRNA)-pos_miRNA_3 << " -- "
										<< (int)strlen(sequence_miRNA)-pos_miRNA_2 << endl;
     cout << "positions with dangle(ncRNA): " << (int)strlen(sequence_miRNA) - pos_miRNA_4_dangle  << " -- " 
          << (int)strlen(sequence_miRNA) - pos_miRNA_1_dangle  << endl;
     //cout << "positions(miRNA ): " << pos_miRNA_1+1  << " -- " << pos_miRNA_4+1  << endl;

     double round_ED_target= (fabs(ED_target[pos_target_1_dangle][pos_target_4_dangle])<DOUBLE_DIFF?0:ED_target[pos_target_1_dangle][pos_target_4_dangle]);
     double round_ED_miRNA = (fabs(ED_miRNA[pos_miRNA_1_dangle][pos_miRNA_4_dangle])<DOUBLE_DIFF?0:ED_miRNA[pos_miRNA_1_dangle][pos_miRNA_4_dangle]);
     cout << "ED target need: " << round_ED_target                  << " kcal/mol" << endl;
     cout << "ED ncRNA  need: " << round_ED_miRNA                   << " kcal/mol" << endl;
     cout << "hybrid energy : " << v-round_ED_target-round_ED_miRNA << " kcal/mol" << endl;
     cout                                                                          << endl;
    }

  cout << "energy: "  << v << " kcal/mol" << endl << endl;

  clear_hydrization_seq();

  return;
}


/*********************************************************************************
* doing all calculations                                                         *
*********************************************************************************/

void calculation()
{
   //init String for suboptimals solutions in every run
   hybrid_pos_target = string((int)strlen(sequence_target),'0');
   pair<bool, double> subopt_hybrid;

   //Computation and Traceback of optimal hybridization   
   //****************************************************
   allocate_global_matrices();  /* space allocation */
   init_global_matrices()    ;  /* init */

   double v=compute_hybrid();
   bool   hybridization_found=traceback(v);

   //Output of optimal hybridization   
   //****************************************************

   if (detailed_output && target_counter==0 && miRNA_counter==0)
       print_input();
   
   if (!hybridization_found || v>=threshold)
     {
       if ((target_file == false) && (miRNA_file == false))
         {
         if (is_ambiguous_sequence(sequence_target) ||
              is_ambiguous_sequence(sequence_miRNA))
            cout << "\n# WARNING: Ambiguous nucleotide(s) 'N' are excluded from all base pairings of the interaction.\n\n";

            cout << "no significant hybridization found" << endl << endl;
         }
       clear_hydrization_seq();
     }
   else
     {
       if ((target_counter > 0) || (miRNA_counter > 0))
	 cout << "=========================" << endl << endl;;

       // Name of target
       cout << target_name << endl;
 
       // Sequence of target
       if (detailed_output)
	 cout << sequence_target << endl;

       // Name of miRNA
       cout << miRNA_name  << endl;

       // Sequence of miRNA
       if (detailed_output) 
         {
	   change_direction(sequence_miRNA);
	   cout << sequence_miRNA << endl;
	   change_direction(sequence_miRNA);
         }

       cout << endl;

         if (is_ambiguous_sequence(sequence_target) ||
              is_ambiguous_sequence(sequence_miRNA))
            cout << "# WARNING: Ambiguous nucleotide(s) 'N' are excluded from all base pairings of the interaction.\n\n";

       print_hybridization(v);

       //Traceback of suboptimal hybridizations and Output
       //****************************************************

       // Find and Output max. n suboptimal hybridizations if below value 'threshold'
       for (int n=0;hybridization_found && n < max_subopt && v < threshold;n++)
         {
           // find next suboptimal hybridization
           subopt_hybrid = find_next_subopt();
           hybridization_found = subopt_hybrid.first;
           v = subopt_hybrid.second;

           if (hybridization_found && v < threshold)
             {
               print_hybridization(v);
             }
           else
            {
              if (n == 0)
                {
                  if ((target_file == false) && (miRNA_file == false))
                    {
                      cout << "no significant suboptimal hybridization found" << endl << endl;
                    }
                 }
              else
               {
                 if ((target_file == false) && (miRNA_file == false))
                   {
                     cout << "no more significant suboptimal hybridizations found" << endl << endl;
                   }
               }
              clear_hydrization_seq();
            }
         }
     }
   
   // print_matrices();  
	
   free_global_matrices();
	
   return;
}

