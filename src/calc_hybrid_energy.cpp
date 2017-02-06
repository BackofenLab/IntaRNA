/***************************************************************************
                       calc_hybrid_energy.cpp  -  main
                             -------------------
    begin                : April-2007
    copyright            : (C) 2007 by Anke Busch & Sven Siebert
    email                : anke.busch@informatik.uni-freiburg.de
 ***************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "basics.h"

using namespace std;


/**********************************************************************************
*    Constants                                                                   *
*          defined in "basics.h"
**********************************************************************************/

const int MAX_BULGE_SIZE = 16;
const int MAX_INTERIOR_ONE_SIDE_SIZE = 16;

const double K_0 = 273.15;         // temperature of 0°C in [K]
const double T_MEASURE  = 310.15;  // temperature in [K] used to measure free energies
const double R  = 1.98717;         // gas constant in [cal/K/mol], value used by Vienna Package 1.7

// const double T                    = 310.15  ; // temperature in Kelvin(37C);K=C+273.15
// const double RT                   = 0.616   ; // Boltzmann gas constant in kcal/mol
// RT=kT=1.380622*10^(-23)Joule/K *310.15K (37C=273.15K+37K=310.15K)


const double MIN_DOUBLE = -1000000.0;
const double MAX_DOUBLE =  1000000.0;
const double DOUBLE_DIFF = 0.00001;

const int MIN_INT = -1000000;
const int MAX_INT = 1000000;

/**********************************************************************************
*                               global variabels                                  *
**********************************************************************************/
double   RT                  ; // product of gas constant and temperature
double   duplex_init         ; // intermolecular initiation free energy
double   terminal_AU         ; // terminal AU penalty (actually non-GC penalty)
double   loop_extrapolation  ; // extrapolation term for large loops (internal, bulge or hairpin) > 30 nt:
                               // dS(T)=dS(30)+loop_extrapolation*ln(n/30)
double   max_ninio_correction; // asymmetric internal loops maximum correction term
double   ninio_correction    ; // asymmetric internal loops Ninio correction term

char*    sequence_target     ;  //input sequence 1
char*    sequence_miRNA      ;  //input sequence 2
int*     int_sequence_target ;  //input sequence 1 as integers (A=0,C=1,G=2,U=3)
int*     int_sequence_miRNA  ;  //input sequence 2 as integers (A=0,C=1,G=2,U=3)

vector< string > targets     ;  //used if sequences are put in via a .fasta file
vector< string > miRNAs      ;  //used if sequences are put in via a .fasta file
int**    che1                ;  //C-hybridization-end in the mRNA
int**    che2                ;  //C-hybridization-end in the sRNA
int**    ches1               ;  //Cseed-hybridization-end in the mRNA
int**    ches2               ;  //Cseed-hybridization-end in the sRNA
double** C                   ;  //C-matrix (~RNAhybrid, includes mfe if the last bases pair)
double** Cseed               ;  //C^seed-matrix (RNAhyrid+seed, includes mfe if the last bases pair)
double** U                   ;  //U-matrix (includes mfe if last bases in one or both sequences are unpaired)
double** ED_target           ;  //matrix of the energy differences of the mRNA (ED = E^unpaired-E^all)
bool     ED_target_computed  ;  //is true if ED_target is computed, which is done after the first valid seed is found
double   ED_target_scaling   ;  //scaling factor for target ED
double** ED_miRNA            ;  //matrix of the energy differences of the sRNA (ED = E^unpaired-E^all)
bool     ED_miRNA_computed   ;  //is true if ED_miRNA is computed, which is done after the first valid seed is found
double   ED_miRNA_scaling    ;  //scaling factor for miRNA ED
double** seed                ;  //matrix, that stores the best seed-hybridizing energy, (seed(i,j) gives the mfe of hybridization starting at (i,j))
int**    seed_unpaired_target;  //matrix, that stores the number of unbound bases in the target that give the mfe in seed
int**    seed_unpaired_miRNA ;  //matrix, that stores the number of unbound bases in the miRNA that give the mfe in seed
int      num_bp_seed         ;  //number of base pairs that the user wants to be in the seed region
int      max_unpaired_target ;  //max. number of unpaired bases in the seed of the target
int      max_unpaired_miRNA  ;  //max. number of unpaired bases in the seed of the target
int      max_unpaired        ;  //max. number of unpaired bases in the seed (sum of unpaired bases in both sequences)
// next line is quick hack, please remove if not needed anymore
double   pu_comb_threshold   ;  // lower end of boxplot for L100, W200 and u=5: 0.0068; u=6: 0.0018; u=7: 0.00045
int      seed_region_start_ncRNA_5;  //start of seed region in ncRNA as passed from 5' end counting from 1
int      seed_region_end_ncRNA_5;  //end of seed region in ncRNA as passed from 5' end counting from 1
int      seed_region_start_ncRNA_3;  //internal start of seed region in ncRNA from 3' end counting from 0
int      seed_region_end_ncRNA_3;  //internal end of seed region in ncRNA given from 3' end counting from 0
bool     seed_region_ncRNA   ;  //is true if seed region in ncRNA is given
bool     window              ;  //is true if sliding window approach is used 
int      window_size         ;  //size of the sliding window
double   threshold           ;  //threshold for output of results below value 'threshold'
bool     max_length          ;  //is true if max. size of unpaired region is specified
int      max_unpaired_length ;  //max. size of the unpaired region during the calculation of ED_...
bool     max_dist            ;  //is true if max. base pair span is specified
int      max_pair_dist       ;  //maximal  base  pair  span used for ED-calculation with RNAplfold
int	 max_subopt	     ;	//max. number of calculated suboptimal hybridizations
int      best_che1           ;  //best C-hybridization-end in the mRNA
int      best_che2           ;  //best C-hybridization-end in the sRNA

bool     detailed_output     ;  //is true, if a detailed output is wanted
bool     use_RNAplfold       ;  //is true, if RNAplfold-procedures are used to compute ED-values in the target
bool     use_RNAup           ;  //is true, if RNAup-procedures are used to compute ED-values in the sRNA
bool     target_file         ;  //true if target file is given 
bool     miRNA_file          ;  //true if miRNA file is given
bool     heuristic           ;  //true if heuristic for hybridization end is used

double***** hybrid_matrix    ;  //stores the hybridizing energies

string   target_name         ;   //name of current target
string   miRNA_name          ;   //name of current miRNA
string   hybrid_sequence_target; //target sequence after hybridization
string   hybrid_sequence_miRNA ; //miRNA sequence after hybridization 
string   bulges_sequence_target; //bulges in target sequence after hybridization
string   bulges_sequence_miRNA ; //bulges in miRNA sequence after hybridization
string   hybrid_pos_target   ;   //positions in target not considered yet for suboptimal hybridizations

int      target_counter      ;  //counts the analyzed targets
int      miRNA_counter       ;  //counts the analyzed miRNAs


/**********************************************************************************
*           stable regions      very stable region(seed)       stable region
*
*                             pos2                       pos3
*          pos1               |                           |                 pos4
*           |                 *****************************                  |
*           ******************************************************************
* ACUAUCGUGCUAGUGCAUGCUGACUGAUGCUGAUCUGAUGCUGAUGCAUCUAUCUAUCUACUAUCUGAUCGAUCGAUCCUAGCUGAGCUAC   
*          ********************************************************************
*          |                                                                  |
*     pos1_dangle                                                        pos4_dangle
*       
**********************************************************************************/

int      pos_target_1  ;
int      pos_target_2  ;
int      pos_target_3  ;
int      pos_target_4  ;
int      pos_target_1_dangle ;
int      pos_target_4_dangle ;
int      pos_miRNA_1   ;
int      pos_miRNA_2   ;
int      pos_miRNA_3   ;
int      pos_miRNA_4   ;
int      pos_miRNA_1_dangle  ;
int      pos_miRNA_4_dangle  ;



/**********************************************************************************
*                        help for calling the program                             *
**********************************************************************************/

void usage()
{
  cout << endl;
  cout << "Synopsis    : IntaRNA [-t fasta_file] [-m fasta_file] [-v energy] [-o]"            << endl;
  cout << "                      [-s number] [-n]  [-h] [-u[1|2] number] [-u number]"         << endl;
  cout << "                      [-p number] [-f pos,pos] [-T temp] [-U] [-P][-w size]"       << endl;
  cout << "                      [-L distance] [-l length] [-a weight] [-b weight]"           << endl;
  cout << "                      [-c threshold] target-RNA-seq binding-RNA-seq"               << endl << endl;
  cout << "Description : finds optimal hybridization location between RNA sequences"          << endl;
  cout << "              target-RNA-seq and binding-RNA-seq according to the sum of 3 energy" << endl;
  cout << "              contributions:"                                                      << endl;
  cout << "                1) energy to open the binding site in the target RNA"              << endl;
  cout << "                2) energy to open the binding site in the binding RNA"             << endl;
  cout << "                3) hybridization energy between both RNAs"                         << endl << endl;;
  cout << "Options     : "                                                                    << endl;
  cout << " === General parameters === "                                                      << endl;
  cout << "              -t fasta_file    : use fasta file of target sequences"               << endl;
  cout << "              -m fasta_file    : use fasta file of binding sequences"              << endl;
  cout << "              -v energy        : outputs all results below energy in kcal/mol"     << endl;
  cout << "              -o               : detailed output"                                  << endl;
  cout << "              -s number        : max. number of calculated suboptimal results"     << endl;
  cout << "                                 (default:0)"                                      << endl;
  cout << "              -n               : use no heuristic for hybridization end"           << endl;
  cout << "                                 (complete approach, more time-consuming)"         << endl;
  cout << "                                 Does not support -s option"                       << endl;
  cout << "              -h               : this help"                                        << endl << endl;
  cout << " === Seed parameters === "                                                         << endl;
  cout << "              -p number        : exact number of paired bases in the seed region"  << endl;
  cout << "                                 (default:6)"                                      << endl;
  cout << "              -u[1|2] number   : max. number of unpaired bases in the seed"        << endl;
  cout << "                                 region in"                                        << endl;
  cout << "                                 1: the first  sequence (default:0)"               << endl;
  cout << "                                 2: the second sequence (default:0)"               << endl;
  cout << "              -u number        : max. number of unpaired bases in the seed"        << endl;
  cout << "                                 region in both sequences (default:0)"             << endl;
  cout << "              -f number,number : search for seed in binding RNA (e.g. ncRNA)"      << endl;
  cout << "                                 in region between positions start,end"            << endl;
  cout << "                                 (given in 5' to 3' direction counting from 1)"    << endl << endl;
  cout << " === RNA folding parameters === "                                                  << endl;
  cout << "              -T temp          : temperature in Celsius (default: 37°C)"           << endl;
  cout << "              -U               : use RNAup to compute ED values of binding RNA"    << endl;
  cout << "                                 (default)"                                        << endl;
  cout << "              -P               : use RNAplfold to compute ED values of target RNA" << endl;
  cout << "                                 (default)"                                        << endl;
  cout << "              -w size          : window size for computation of ED values"         << endl;
  cout << "                                 with RNAplfold (default: length of target RNA)"   << endl;
  cout << "              -L distance      : max. distance of two paired bases for"            << endl;
  cout << "                                 computation of ED values with RNAplfold"          << endl;
  cout << "                                 (default: window size)"                           << endl;
  cout << "              -l length        : max. length of hybridized region, mainly used "   << endl;
  cout << "                                 for efficient computation (default: window size)" << endl;
  cout << "              -a weight        : weight for ED values of target RNA in energy     "<< endl;
  cout << "                                 (default: 1.0)"                                   << endl;
  cout << "              -b weight        : weight for ED values of binding RNA in energy    "<< endl;
  cout << "                                 (default: 1.0)"                                   << endl;
  cout << "              -c threshold     : threshold for seed accessibility, requires u=0   "<< endl;
  cout << "                                 EXPERIMENTAL FEATURE, (default: -1.0)"            << endl << endl;
  cout << "Example     :-find optimal hybridization location between two RNAs:"               << endl;
  cout << "              IntaRNA CCCCCCCCCCCCGGG GGGGGGCCC"                                   << endl;
  cout << "              Output: "                                                            << endl;
  cout << "                        5'-CCCCCC      GGG-3'"                                     << endl;
  cout << "                                 CCCCCC"                                           << endl;
  cout << "                                 GGGGGG"                                           << endl;
  cout << "                           3'-CCC      -5'"                                        << endl;
  cout << "                      energy: -11.4152  kcal/mol"                                  << endl << endl;
  cout << "Contact     : intarna@informatik.uni-freiburg.de"                                  << endl;
  cout << endl;
}

/*
void usage_help(char* name)
{
   cout << "\ncall: " << name << " \"sequence1\" \"sequence2\"\n";
   cout << "             [-p number of base paired in the seed region, default: 4]\n";
   cout << "             [-B max. number of unpaired bases in the seed region of seq. 1, default: 0]\n";
   cout << "             [-b max. number of unpaired bases in the seed region of seq. 2, default: 0]\n";
   cout << "             [-U max. number of unpaired bases in the seed region of both seq's, default: 0] \n";
   cout << "             [-w window size for sliding window approach, default: sliding window appr. not used]\n";
   cout << "             [-l max. length of an unpaired region (that is considered)]\n";
   cout << "             [-T temperature in Celsius, default: 37 Celsius]\n";
   cout << "             [-o] if a detailed output is wanted\n";
   cout << "\n";
   cout << endl;
   exit(1);
}
*/



/**********************************************************************************
*                                 Main-Function                                   *
**********************************************************************************/

int main(int argc, char *argv[])
{
   do_backtrack        = 1; // must be set for RNAup
   fold_constrained    = 1;
   
   /*default values*/

// energy parameter from mfold Version 3.0 (Mathews,Sabina, Zuker & Turner:JMB1999)- do not change!
   RT                  = R*T_MEASURE/1000.0; // in [kcal/mol]
   duplex_init         = 4.10;
   terminal_AU         = 0.50;
   loop_extrapolation  = 1.079;
   max_ninio_correction = 3.00;
   ninio_correction    = 0.50;
// ************************************************************************************************

   ED_target_computed  = false;
   ED_miRNA_computed   = false;

   ED_target_scaling   = 1.0;
   ED_miRNA_scaling    = 1.0;

   num_bp_seed         = 6;
   max_unpaired_target = 0;
   max_unpaired_miRNA  = 0;
   max_unpaired        = 0;
   pu_comb_threshold   = -1.0;
   seed_region_start_ncRNA_5  = -1;
   seed_region_start_ncRNA_3  = -1;
   seed_region_end_ncRNA_5  = -1;
   seed_region_end_ncRNA_3  = -1;
   seed_region_ncRNA   = false;
   window              = false;
   window_size         = -1;
   max_length          = false;
   max_unpaired_length = -1;
   max_dist            = false;
   max_pair_dist       = -1;
   detailed_output     = false;
   threshold           = MAX_DOUBLE;
   max_subopt	       = 0;
   best_che1           = -1;
   best_che2           = -1;
   target_name         = "";
   miRNA_name          = "";
   use_RNAplfold       = true;
   use_RNAup           = true;
   heuristic           = true;
   
   target_file = false;
   char* target_file_name = NULL;
   
   miRNA_file = false;
   char* miRNA_file_name = NULL;
   
   target_counter = 0;
   miRNA_counter  = 0;

   /**********************************************
   *  reading the options                        *
   **********************************************/   
     
   int opt;
   int option_index = 0;
   
   static struct option long_options[] =
   {
      /* These options are distinguished by their indices. */
      {"u1",  required_argument, 0, 'x'},
      {"u2",  required_argument, 0, 'y'},
	  {0,0,0,0}
   };
   
   while ((opt = getopt_long_only(argc, argv, "x:y:u:p:f:T:w:l:L:a:b:c:v:s:t:m:oPnhU",long_options, &option_index)) != -1)
   {
      switch (opt)
	  {
	
	     case 'x':  max_unpaired_target = atoi(optarg);
                    //cout << "max_unpaired_target: " << max_unpaired_target << endl;
                    break;
					
         case 'y':  max_unpaired_miRNA = atoi(optarg);
                    //cout << "max_unpaired_miRNA : " << max_unpaired_miRNA << endl;
                    break;
					
		 case 'u':  max_unpaired = atoi(optarg);
                    //cout << "max_unpaired       : " << max_unpaired << endl;
                    break;

 	 /*case 'u': int seqnum      =atoi(argv[optind]);
	           int num_unpaired=atoi(argv[optind+1]);
 		   optind+=2;
 	           if (seqnum==1)
 		     max_unpaired_target = num_unpaired;
 		   else if (seqnum==2)
 		     max_unpaired_miRNA = num_unpaired;
 		   else if (seqnum==3)
 		     max_unpaired = num_unpaired;
 		   else 
 		    {
 		      cerr << "no valid sequence identifier" << endl;
 		      exit(1);
 		    }
 		  break;*/
         
	 case 'p': num_bp_seed = atoi(optarg);
	           break;
	 case 'f': // passed positions count from 5' with 1, but internal indices count from 3' with 0
			   sscanf(optarg,"%d,%d",&seed_region_start_ncRNA_5,&seed_region_end_ncRNA_5);
               seed_region_ncRNA = true;
	           break;
	 case 'T': temperature = atof(optarg); // in [°C]
                   RT = R*(temperature+K_0)/1000.0; // in [kcal/mol]
                   scale_parameters();
		   break;
	 case 'w': window_size = atoi(optarg);
	           window = true;
		   break;
	 case 'l': max_unpaired_length = atoi(optarg);
	           max_length = true;
		   break;
	 case 'L': max_pair_dist = atoi(optarg);
	           max_dist = true;
		   break;
	 case 'a': ED_target_scaling = atof(optarg);
		   break;
	 case 'b': ED_miRNA_scaling = atof(optarg);
		   break;
	 case 'c': pu_comb_threshold = atof(optarg);
		   break;
	 case 'v': threshold = atoi(optarg);
                   break;
	 case 's': max_subopt= (heuristic?atoi(optarg):0);
                   break;
	 case 't': target_file_name = optarg;
	           target_file = true;
		   break;
	 case 'm': miRNA_file_name = optarg;
	           miRNA_file = true;
		   break;
	 case 'o': detailed_output = true;
	           break;
	 case 'P': use_RNAplfold = true;
	           break;
	 case 'U': use_RNAup = true;
						// must be set to calculate the base pair probability matrix
						 do_backtrack = 1;
	           break;
         case 'n': heuristic = false;
                   max_subopt = 0;
                   break;
	 case 'h': usage();
	           exit(0);
         case '?': if (isprint (optopt) )
                      cerr << "Unknown option `-"          << optopt << "`" << endl;
                   else
                      cerr << "Unknown option character `" << optopt << "`" << endl;
               return 1;
	 default : abort ();
	 }
  }
   
   /*****************************************************************
   * if RNAplfold is used, max_length and window must also be given *
   *****************************************************************/
   
   /*if ((use_RNAplfold == true) && ((max_length == false) || (window == false)))
   {
      cerr << "Window size (-w) and the max. length of the hybridized region (-l) must be given when using RNAplfold (-P)!"   << endl;
      exit(1);
   }*/
   
   /**********************************************
    *  reading the sequences                      *
    **********************************************/
   
   if ((target_file == false) && (miRNA_file == false))
     {
       if (optind < argc)
	 {
	   sequence_target = argv[optind++];
	   if (optind < argc)
	     {
	       sequence_miRNA = argv[optind++];
	       if (optind < argc)
		 {
		   for (int i=optind;i<argc;i++)
		     cerr << "Non-option argument: " << argv[i] << endl;
		   exit(1);
		 }
	     }
	   else
	     {
	       cerr << "The second sequence is missing!\n";
	       exit(1);
	     }
	 }
       else
	 {
	   cerr << "Both RNA sequences are missing!\n";
	   exit(1);
	 }
     }
   
   else if ((target_file == true) && (miRNA_file == false))
     {
       if (optind < argc)
	 {
	   sequence_miRNA = argv[optind++];
	   if (optind < argc)
	     {
	       for (int i=optind;i<argc;i++)
		 cerr << "Non-option argument: " << argv[i] << endl;
	       exit(1);
	     }
	 }
       else
	 {
	   cerr << "The second sequence is missing!\n";
	   exit(1);
	 }
     }
   
   else if ((target_file == false) && (miRNA_file == true))
     {
       if (optind < argc)
	 {
	   sequence_target = argv[optind++];
	   if (optind < argc)
	     {
	       for (int i=optind;i<argc;i++)
		 cerr << "Non-option argument: " << argv[i] << endl;
	       exit(1);
	     }
	 }
       else
	 {
	   cerr << "The first sequence is missing!\n";
	   exit(1);
	 }
     }
   else
     {
       if (optind < argc)
	 {
	   for (int i=optind;i<argc;i++)
	     cerr << "Non-option argument: " << argv[i] << endl;
	   exit(1);
	 }
     }
   
   
   
   // reading in files
   bool okay;   
   if (target_file == true)
     {
       okay = readin_fasta(target_file_name, targets);
       if (!okay)
	 {
	   cerr << "No valid sequence in file " << target_file_name << "!\n";
	   exit(1);
	 }
     }

   if (miRNA_file == true)
     {
       okay = readin_fasta(miRNA_file_name, miRNAs);
       if (!okay)
	 {
	   cerr << "No valid sequence in file " << miRNA_file_name << "!\n";
	   exit(1);
	 }
     }
	
   /***************************************************************************
   * different cases of calculations (depending on the input format)          *
   ***************************************************************************/
	
	
   if ((target_file == false) && (miRNA_file == false))
     {
       //test, whether the sequences are valid and sets it to upper case letters
       if (!check_sequence_and_set2upper(sequence_target))
	 {
	   cerr << "\nSequence 1 is not valid!\n\n";
	   exit(1);
	 }
       if (!check_sequence_and_set2upper(sequence_miRNA))
	 {
	   cerr << "\nSequence 2 is not valid!\n\n";
	   exit(1);
	 }

       target_name = ">target RNA";
       miRNA_name  = ">ncRNA";
       //hybrid_pos_target = string((int)strlen(sequence_target),'0');  //moved to calculation()
       
       //change direction of the second sequence
       change_direction(sequence_miRNA);
       // passed ncRNA positions counts from 5' with 1, but internal indices count from 3' with 0
       if (seed_region_ncRNA) {
          seed_region_start_ncRNA_3 = (int)strlen(sequence_miRNA) - seed_region_end_ncRNA_5;
          seed_region_end_ncRNA_3 = (int)strlen(sequence_miRNA) - seed_region_start_ncRNA_5;
       }

         /***************************************************************************
         * test, whether the values of the input options relevant for              *
         * hybridization are valid                                                 *
         ***************************************************************************/
         
         if (!test_options_hybridization())
         exit(1);	

       init_target();
       init_miRNA();
       
//    cout << "window: " << window << endl;
//    cout << "window_size: " << window_size << endl;
//    cout << "max_length: " << max_length << endl;
//    cout << "max_unpaired_length: " << max_unpaired_length << endl;

       calculation();

       free_target();
       free_miRNA();
     }
   
   else if ((target_file == true) && (miRNA_file == false))
     {
       if (!check_sequence_and_set2upper(sequence_miRNA))
	 {
	   cerr << "\nSequence " << sequence_miRNA << " is not valid!\n\n";
	   exit(1);
	 }
       
       //change direction of the second sequence
       change_direction(sequence_miRNA);
       // passed ncRNA positions counts from 5' with 1, but internal indices count from 3' with 0
       if (seed_region_ncRNA) {
          seed_region_start_ncRNA_3 = (int)strlen(sequence_miRNA) - seed_region_end_ncRNA_5;
          seed_region_end_ncRNA_3 = (int)strlen(sequence_miRNA) - seed_region_start_ncRNA_5;
       }
       miRNA_name  = ">ncRNA";
       init_miRNA();
       for (int seq = 1; seq < (int)targets.size(); seq+=2)
	 {
	   sequence_target = (char*)malloc(sizeof(char)*((int)targets[seq].size()+1));
	   for (int i=0; i<(int)targets[seq].size(); i++)
	     sequence_target[i] = targets[seq][i];
	   
	   sequence_target[(int)targets[seq].size()] = '\0';

	   target_name = targets[seq-1];

           //hybrid_pos_target = string((int)strlen(sequence_target),'0');  //moved to calculation()

            /***************************************************************************
            * test, whether the values of the input options relevant for              *
            * hybridization are valid                                                 *
            ***************************************************************************/
            
            if (!test_options_hybridization())
            exit(1);

           init_target();

	   calculation();
	   target_counter++;

           free_target();
	   free(sequence_target);
	 }
       free_miRNA();
     }
   
   else if ((target_file == false) && (miRNA_file == true))
     {
       if (!check_sequence_and_set2upper(sequence_target))
       {
         cerr << "\nSequence " << sequence_target << " is not valid!\n\n";
         exit(1);
       }
       //hybrid_pos_target = string((int)strlen(sequence_target),'0');    //moved to calculation()

       init_target();
       for (int mi = 1; mi < (int)miRNAs.size(); mi+=2)
	 {
	   sequence_miRNA = (char*)malloc(sizeof(char)*((int)miRNAs[mi].size()+1));
	   for (int i=0; i<(int)miRNAs[mi].size(); i++)
	     sequence_miRNA[i] = miRNAs[mi][i];
	   sequence_miRNA[(int)miRNAs[mi].size()] = '\0';
	   
	   //change direction of the second sequence
	   change_direction(sequence_miRNA);
       // passed ncRNA positions counts from 5' with 1, but internal indices count from 3' with 0
       if (seed_region_ncRNA) {
          seed_region_start_ncRNA_3 = (int)strlen(sequence_miRNA) - seed_region_end_ncRNA_5;
          seed_region_end_ncRNA_3 = (int)strlen(sequence_miRNA) - seed_region_start_ncRNA_5;
       }

	   target_name = ">target RNA";
	   miRNA_name  = miRNAs[mi-1];

            /***************************************************************************
            * test, whether the values of the input options relevant for              *
            * hybridization are valid                                                 *
            ***************************************************************************/
            
            if (!test_options_hybridization())
            exit(1);

           init_miRNA();

	   calculation();
	   miRNA_counter++;
           free_miRNA();
	   free(sequence_miRNA);
	 }
       free_target();
     }
   else 
     {
       for (int seq = 1; seq < (int)targets.size(); seq+=2)
	 {
	   sequence_target = (char*)malloc(sizeof(char)*((int)targets[seq].size()+1));
	   for (int i=0; i<(int)targets[seq].size(); i++)
	     sequence_target[i] = targets[seq][i];
	   sequence_target[(int)targets[seq].size()] = '\0';

	   target_name = targets[seq-1];
	   //hybrid_pos_target = string((int)strlen(sequence_target),'0');  //moved to calculation()

           init_target();
	   for (int mi = 1; mi < (int)miRNAs.size(); mi+=2)
	     {
	       sequence_miRNA = (char*)malloc(sizeof(char)*((int)miRNAs[mi].size()+1));
	       for (int i=0; i<(int)miRNAs[mi].size(); i++)
		 sequence_miRNA[i] = miRNAs[mi][i];
	       sequence_miRNA[(int)miRNAs[mi].size()] = '\0';
	       
	       //change direction of the second sequence
	       change_direction(sequence_miRNA);
           // passed ncRNA positions counts from 5' with 1, but internal indices count from 3' with 0
          if (seed_region_ncRNA) {
              seed_region_start_ncRNA_3 = (int)strlen(sequence_miRNA) - seed_region_end_ncRNA_5;
              seed_region_end_ncRNA_3 = (int)strlen(sequence_miRNA) - seed_region_start_ncRNA_5;
          }

	       miRNA_name  = miRNAs[mi-1];

               /***************************************************************************
               * test, whether the values of the input options relevant for              *
               * hybridization are valid                                                 *
               ***************************************************************************/
               
               if (!test_options_hybridization())
               exit(1);

               init_miRNA();
	       
	       calculation();
	       miRNA_counter++;
               free_miRNA();
	       free(sequence_miRNA);
	     }
	   target_counter++;
           free_target();
	   free(sequence_target);
	   miRNA_counter = 0;
	 }
     }
   
   return 0;
}




