#include "basics.h"

/**********************************************************
 translates an integerBase to a characterBase (nucleotide)
**********************************************************/

char int2char(int base)
{
   if (base == 0) return 'A';
   else if (base == 1) return 'C';
   else if (base == 2) return 'G';
   else if (base == 3) return 'U';
   else if (base == 4) return 'N';
   else
   {
      cerr << "There's something wrong in your int-base!";
      exit(1);
   }
}

/**********************************************************
 translates an integerSeq. to a characterSeq. (nucleotide)
**********************************************************/

char* int2char(int* num_seq, int size)
{
   char* seq = (char*) malloc(sizeof(char)*size);
   for (int i=0; i<size; i++)
   {
      if (num_seq[i] == 0) seq[i]='A';
      else if (num_seq[i] == 1) seq[i]='C';
      else if (num_seq[i] == 2) seq[i]='G';
      else if (num_seq[i] == 3) seq[i]= 'U';
      else if (num_seq[i] == 4) seq[i]= 'N';
      else
      {
         cerr << "There's something wrong in your int-sequence!";
         exit(1);
      }
   }
   return seq;
}


/**********************************************************
 translates a characterBase to an integerBase (nucleotide)
**********************************************************/

int char2int_base(char base)
{
   if (toupper(base) == 'A') return 0;
   else if (toupper(base) == 'C') return 1;
   else if (toupper(base) == 'G') return 2;
   else if ((toupper(base) == 'U') || (toupper(base) == 'T')) return 3;
   else if (toupper(base) == 'N') return 4;
   else
   {
      printf("There's something wrong in your char-base!\n");
      exit(1);
   }
}

/**********************************************************
 translates a characterSeq. to an integerSeq. (nucleotide)
**********************************************************/

int* char2int(char* seq)
{
   int* num_seq = (int*) malloc(sizeof(int)*strlen(seq));
   for (int i=0; i<(int)strlen(seq); i++)
      num_seq[i] = char2int_base(seq[i]);
   return num_seq;
}



/*********************************************************
 changes the direction of the sequence (5'-3' => 3'-5')
*********************************************************/
void change_direction(char* sequence)
{
   int len = (int)strlen(sequence);
   char reverse_sequence[len+1];
   
   for (int i=0; i<len; i++)
      reverse_sequence[len-1-i] = sequence[i];
   reverse_sequence[len] = '\0';

   strncpy(sequence, reverse_sequence, len+1);

   return;
}



/**********************************************************
 translates the two bases (integer) of a BP to an integer
 value for the pair (0,1,2,3,4,5)
**********************************************************/

int BP2int(int b1, int b2)
{
   int kind_of_BP = 0;
   if ( b1==0 && b2==3 ) kind_of_BP = 0;
   else if ( b1==1 && b2==2 ) kind_of_BP = 1;
   else if ( b1==2 && b2==1 ) kind_of_BP = 2;
   else if ( b1==3 && b2==0 ) kind_of_BP = 3;
   else if ( b1==2 && b2==3 ) kind_of_BP = 4;
   else if ( b1==3 && b2==2 ) kind_of_BP = 5;
   else
   {
      cerr << "Wrong Basepair (" << int2char(b1) << " - " << int2char(b2) << ") !" << endl;
      exit(1);
   }
   return kind_of_BP;
}


/**********************************************************
 translates a BP (given as an integer) to two integerBases
**********************************************************/

void BP2_2(int bp, int &base_i, int &base_j)
{
   if (bp == 0) {base_i=0; base_j=3;}
   else if (bp == 1) {base_i=1; base_j=2;}
   else if (bp == 2) {base_i=2; base_j=1;}
   else if (bp == 3) {base_i=3; base_j=0;}
   else if (bp == 4) {base_i=2; base_j=3;}
   else if (bp == 5) {base_i=3; base_j=2;}
   else
   {
      cerr << "Wrong BasepairNr.!" << endl;
      exit(1);
   }
}


/**********************************************************
checks whether 2 bases can pair
**********************************************************/

bool PairTest(int b1, int b2)
{
   if (( b1==0 && b2==3 ) || ( b1==1 && b2==2 ) || ( b1==2 && b2==1 ) || ( b1==3 && b2==0 ) || ( b1==2 && b2==3 ) || ( b1==3 && b2==2 ))
      return true;
   else
      return false;
}

/**********************************************************
checks whether 2 bases can pair
**********************************************************/

bool PairTest(char b1, char b2)
{
   if (( b1=='A' && b2=='U' ) || ( b1=='C' && b2=='G' ) || ( b1=='G' && b2=='C' ) || ( b1=='U' && b2=='A' ) || ( b1=='G' && b2=='U' ) || ( b1=='U' && b2=='G' ))
      return true;
   else
      return false;
}

/*******************************************************************
 tests whether the sequence consists of valid nucleotides (A,C,G,T/U)
 allow also nucleotide N (but should be forbidden for base pairings)
*******************************************************************/

bool check_sequence_and_set2upper(char* seq)
{
   for (int i=0; i<(int)strlen(seq); i++)
   {
      seq[i] = toupper(seq[i]);
      if (seq[i] == 'T')
         seq[i] = 'U';
      if ((seq[i] != 'A') && (seq[i] != 'C') && (seq[i] != 'G') && (seq[i] != 'U') && (seq[i] != 'N'))
         return false;
   }
   return true;
}

bool check_sequence_and_set2upper(string& seq)
{
   for (int i=0; i<(int)seq.size(); i++)
   {
      seq[i] = toupper(seq[i]);
      if (seq[i] == 'T')
         seq[i] = 'U';
      if ((seq[i] != 'A') && (seq[i] != 'C') && (seq[i] != 'G') && (seq[i] != 'U') && (seq[i] != 'N'))
         return false;
   }
   return true;
}

/*******************************************************************
 tests whether the sequence contains nucleotides N
*******************************************************************/

bool is_ambiguous_sequence(char* seq)
{
   for (int i=0; i<(int)strlen(seq); i++)
   {
      if (seq[i] == 'N')
         return true;
   }
   return false;
}

/**********************************************************
since addition and substraction of two float or double 
values is often not exact enough, calculations using float
or double values were rounded to 3 positions after decimal
point: the number is multiplied by 1000, 0.5 is added (if
the number is negative, 0.5 is substracted), then turn the
number to an integer, i.e. the remaining positions after
decimal point are cut. 
(The procedure is similar to rounding but 0.5 is added
before.)
**********************************************************/

double Double_Round(double x)
{
   int x_int;
   if (x>=0)
      x_int = (int)((x*1000)+0.5);
   else
      x_int = (int)((x*1000)-0.5);

   return (double) x_int/1000;
}



/************************************************************
 finds the minimum of two integers
************************************************************/

int Minimum(int a, int b)
{
   if (a<b)
      return a;
   else
      return b;
}

/************************************************************
 finds the minimum of two double values
************************************************************/

double Minimum(double a, double b)
{
   if (a<b)
      return a;
   else
      return b;
}


/************************************************************
 finds the maximum of two integers
************************************************************/

int Maximum(int a, int b)
{
   if (a>b)
      return a;
   else
      return b;
}



/************************************************************
 returns the sum of a and b, if either a or b is MAX_DOUBLE,
 it returns MAX_DOUBLE
************************************************************/

double Sum_or_MAX_DOUBLE(double a, double b)
{
   if ((a==MAX_DOUBLE) || (b==MAX_DOUBLE))
      return MAX_DOUBLE;
   else
      return a+b;
}



/************************************************************
 sums up all values of a vector
************************************************************/

int SumVec(int* vec, int size)
{
   int sum = 0;
   for (int i=0; i<size; i++)
      sum += vec[i];
   return sum;
}

/************************************************************
 finds the minimal value in a vector (returns the value and
 its coordinate)
************************************************************/

double* MiniVec(double* vec, int size)
{
   double* min = (double*) malloc(sizeof(double)*2);
   min[0] = MAX_DOUBLE; //coord.
   min[1] = MAX_DOUBLE; //value
   for ( int i=0; i<size; i++)
      if (vec[i] < min[1])
      {
         min[1] = vec[i];
         min[0] = i;
      }
   return min;
}

/************************************************************
 finds the maximal value in a vector (returns the value and
 its coordinate)
************************************************************/

double* MaxiVec(double* vec, int size)
{
   double* max = (double*) malloc(sizeof(double)*2);
   max[0] = MIN_INT; //coord.
   max[1] = MIN_INT; //value
   for (int i=0; i<size; i++)
      if (vec[i] > max[1])
      {
         max[1] = vec[i];
         max[0] = i;
      }
   return max;
}



