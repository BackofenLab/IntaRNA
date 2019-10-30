
#ifndef INTARNA_INTERACTIONENERGYVIENNA_H_
#define INTARNA_INTERACTIONENERGYVIENNA_H_

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/VrnaHandler.h"

extern "C" {
	#include <ViennaRNA/utils.h>
	#include <ViennaRNA/fold_vars.h>
	#include <ViennaRNA/model.h>
	#include <ViennaRNA/params.h>
	#include <ViennaRNA/loop_energies.h>
}
#ifndef VIENNA_RNA_PAIR_MAT_H
#define VIENNA_RNA_PAIR_MAT_H
extern "C" {
	#include <ViennaRNA/pair_mat.h>
}
#endif

#include <boost/numeric/ublas/triangular.hpp>

#define Evrna_2_E( e ) ( static_cast<E_type>(e) )

namespace IntaRNA {

// http://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/RNAlib-2.3.0.pdf

/**
 * Implements an energy interface based on free energy estimates computed
 * with the Vienna RNA package.
 *
 * @author Martin Mann 2014
 */
class InteractionEnergyVrna: public InteractionEnergy {

public:



	/**
	 * Construct energy utility object given the accessibility ED values for
	 * two sequences.
	 *
	 * @param accS1 accessibility of the first sequence
	 * @param accS2 accessibility of the second sequence
	 * @param vrnaHandler the VRNA parameter handler to be used
	 * @param maxInternalLoopSize1 maximal number of enclosed unpaired positions
	 *          between two intermolecular base pairs in sequence 1, ie it holds
	 *          for an intermolecular loop closed by base pairs (i1,i2) and
	 *          (j1,j2) : (j1-i1+1) <= maxInternalLoopSize
	 * @param maxInternalLoopSize2 maximal number of enclosed unpaired positions
	 *          between two intermolecular base pairs in sequence 2, ie it holds
	 *          for an intermolecular loop closed by base pairs (i1,i2) and
	 *          (j1,j2) : (j2-i2+1) <= maxInternalLoopSize
	 * @param initES whether or not ES values are to be computed
	 * @param energyAdd when computing the overall energy via getE(), this term
	 *          is always added; thus it defines a shift of the energy spectrum
	 *          as e.g. needed when computing predictions with accessibility
	 *          constraints
	 * @param energyWithDangles whether or not dangling end contributions are
	 *          considered within overall energies
	 * @param internalLoopGU whether or not GU base pairs are allowed within
	 *          internal loops
	 */
	InteractionEnergyVrna( const Accessibility & accS1
					, const ReverseAccessibility & accS2
					, VrnaHandler &vrnaHandler
					, const size_t maxInternalLoopSize1 = 16
					, const size_t maxInternalLoopSize2 = 16
					, const bool initES = false
					, const E_type energyAdd = Ekcal_2_E(0.0)
					, const bool energyWithDangles = true
					, const bool internalLoopGU = true
				);

	virtual ~InteractionEnergyVrna();


	/**
	 * Provides the ensemble energy (ES) of all intramolecular substructures
	 * that can be formed within a given region of sequence 1 under the
	 * assumption that the region is part of an (intermolecular) multiloop,
	 * i.e. at least one base pair is formed by each substructure.
	 *
	 * If no structure can be formed within the region, E_INF is returned.
	 *
	 * @param i1 the start of the structured region of seq1
	 * @param j1 the end of the structured region of seq1
	 * @return the ES value for [i1,j1] or E_INF if no intramolecular
	 *         structure can be formed
	 */
	virtual
	E_type
	getES1( const size_t i1, const size_t j1 ) const;

	/**
	 * Provides the ensemble energy (ES) of all intramolecular substructures
	 * that can be formed within a given region of sequence 2 under the
	 * assumption that the region is part of an (intermolecular) multiloop,
	 * i.e. at least one base pair is formed by each substructure.
	 *
	 * If no structure can be formed within the region, E_INF is returned.
	 *
	 * @param i2 the start of the structured region of seq2
	 * @param j2 the end of the structured region of seq2
	 * @return the ES value for [i2,j2] or E_INF if no intramolecular
	 *         structure can be formed
	 */
	virtual
	E_type
	getES2( const size_t i2, const size_t j2 ) const;

	/**
	 * Provides the energy contribution for a given number of unpaired
	 * nucleotides under the
	 * assumption that the region is part of an (intermolecular) multiloop.
	 *
	 * @param numUnpaired the number of unpaired bases
	 * @return the energy contribution of the given number of unpaired bases
	 *         within an intramolecular multiloop
	 */
	virtual
	E_type
	getE_multiUnpaired( const size_t numUnpaired ) const;

	/**
	 * Provides the energy contribution/penalty of the helix repesented by the
	 * interaction right of a multi-site gap starting with base pair (j1,j2)
	 *
	 * @param j1 the end of the gap in seq1, ie the first base paired in the
	 *           interaction site to the right of the gap
	 * @param j2 the end of the gap in seq2, ie the first base paired in the
	 *           interaction site to the right of the gap
	 *
	 * @return the energy contribution/penalty of the intermolecular helix
	 *         within an intramolecular multiloop
	 */
	virtual
	E_type
	getE_multiHelix( const size_t j1, const size_t j2 ) const;

	/**
	 * Provides the energy contribution/penalty for closing an intermolecular
	 * multiloop on the left of a multi-site gap.
	 *
	 * @return the energy contribution/penalty of the intermolecular helix
	 *         within an intramolecular multiloop
	 */
	virtual
	E_type
	getE_multiClosing() const;

	/**
	 * Provides the duplex initiation energy.
	 *
	 * @return the energy for duplex initiation
	 */
	virtual
	E_type
	getE_init() const;

	/**
	 * Computes the energy estimate for the interaction loop region closed by
	 * the intermolecular base pairs (i1,i2) and (j1,j2) where the regions
	 * [i1,j1] and [i2,j2] are considered unpaired.
	 * The energy estimate is derived via the Vienna RNA package loop energies
	 * or is E_INF if the internal loop size exceeds
	 * the allowed maximum (see constructor).
	 *
	 * @param i1 the index of the first sequence (<j1) interacting with i2
	 * @param j1 the index of the first sequence (>i1) interacting with j2
	 * @param i2 the index of the second sequence (<j2) interacting with i1
	 * @param j2 the index of the second sequence (>i2) interacting with j1
	 *
	 * @return energy in kcal/mol for the loop closed by (i1,i2)
	 *         or
	 *         E_INF if the allowed loop size is exceeded or no valid internal loop boundaries
	 */
	virtual
	E_type
	getE_interLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;


	/**
	 * Computes the dangling end energy penalties for the left side
	 * (i1-1 and i2-1) of the interaction closed by the intermolecular
	 * base pair (i1,i2).
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return the dangling end penalty for the left side of the interaction
	 */
	virtual
	E_type
	getE_danglingLeft( const size_t i1, const size_t i2 ) const;


	/**
	 * Computes the dangling end energy penalties for the right side
	 * (j1+1 and j2+1) of the interaction closed by the intermolecular
	 * base pair (j1,j2).
	 *
	 * @param j1 the index of the first sequence interacting with j2
	 * @param j2 the index of the second sequence interacting with j1
	 *
	 * @return the dangling end penalty for the right side of the interaction
	 */
	virtual
	E_type
	getE_danglingRight( const size_t j1, const size_t j2 ) const;


	/**
	 * Provides the penalty for closing an interaction with the given
	 * base pair on the "left side" (i1 = 5' end of seq1 of the interaction)
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return the loop closure penalty for the left side of the interaction
	 */
	virtual
	E_type
	getE_endLeft( const size_t i1, const size_t i2 ) const;

	/**
	 * Provides the penalty for closing an interaction with the given
	 * base pair on the "right side" (j1 = 3' end of seq1 of the interaction)
	 *
	 * @param j1 the index of the first sequence interacting with j2
	 * @param j2 the index of the second sequence interacting with j1
	 *
	 * @return the loop closure penalty for the right side of the interaction
	 */
	virtual
	E_type
	getE_endRight( const size_t j1, const size_t j2 ) const;

	/**
	 * Provides the overall ensemble energy for sequence 1
	 * given its accessibility constraints
	 * @return Eall(constraint-conform intra-molecular structures for seq1)
	 */
	virtual
	E_type
	getEall1() const;

	/**
	 * Provides the overall ensemble energy for sequence 2
	 * given its accessibility constraints
	 * @return Eall(constraint-conform intra-molecular structures for seq2)
	 */
	virtual
	E_type
	getEall2() const;

	/**
	 * Provides the overall partition function for sequence 2
	 * given its accessibility constraints
	 * @return Z(constraint-conform intra-molecular structures for seq2)
	 */
	virtual
	Z_type
	getRT() const;

protected:


	//! Vienna RNA package : folding model to be used for the energy computation
	vrna_md_t foldModel;

	//! Vienna RNA package : folding parameters to be used for the energy
	//! computation
	vrna_param_t * foldParams;

	//! the RT constant to be used for Boltzmann weight computations
	Z_type RT;

	//! base pair code for (C,G)
	const int bpCG;

	//! base pair code for (G,C)
	const int bpGC;

	//! matrix to store ES values (upper triangular matrix)
	typedef boost::numeric::ublas::triangular_matrix<E_type, boost::numeric::ublas::upper> EsMatrix;

	//! the ES values for seq1 if computed (otherwise NULL)
	EsMatrix * esValues1;

	//! the ES values for seq2 if computed (otherwise NULL)
	EsMatrix * esValues2;

	//! ensemble energy of intra-molecular structures of seq1
	mutable E_type Eall1;

	//! ensemble energy of intra-molecular structures of seq2
	mutable E_type Eall2;

	/**
	 * Checks whether or not a given base pair is a GC base pair
	 * @param i1 the index in the first sequence
	 * @param i2 the index in the second sequence
	 * @return true if (seq1(i1),seq2(i2)) is (G,C) or (C,G)
	 */
	bool
	isGC( const size_t i1, const size_t i2 ) const;

	/**
	 * Computes the ES values and fills esValues container
	 * @param acc the accessibility object for the sequence to compute the ES values for
	 * @param esToFill the container to write the ES values to
	 */
	void
	computeES( const Accessibility & acc, EsMatrix & esToFill );

	/**
	 * Computes the ensemble energy of all intra-molecular structures that
	 * are conform to the accessibility constraints
	 * @param acc the Accessibility object to access sequence and constraints
	 * @return the computed ensemble energy
	 */
	E_type
	computeIntraEall( const Accessibility & acc ) const;
};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getE_init() const
{
	// init term is sequence independent
	return Evrna_2_E(foldParams->DuplexInit);
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getE_endLeft( const size_t i1, const size_t i2 ) const
{
	// VRNA non-GC penalty
	return Evrna_2_E( isGC(i1,i2) ? 0.0 : foldParams->TerminalAU );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getE_endRight( const size_t j1, const size_t j2 ) const
{
	// VRNA non-GC penalty
	return Evrna_2_E( isGC(j1,j2) ? 0.0 : foldParams->TerminalAU );
}

////////////////////////////////////////////////////////////////////////////

inline
Z_type
InteractionEnergyVrna::
getRT() const
{
	return RT;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getE_interLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	// if valid internal loop
	if ( isValidInternalLoop(i1,j1,i2,j2) ) {
		assert( i1!=j1 && i2!=j2 );
		// Vienna RNA : compute internal loop / stacking energy for base pair [i1,i2]
		return Evrna_2_E(E_IntLoop(	(int)j1-i1-1	// unpaired region 1
							, (int)j2-i2-1	// unpaired region 2
							, BP_pair[accS1.getSequence().asCodes().at(i1)][accS2.getSequence().asCodes().at(i2)]	// type BP (i1,i2)
							, BP_pair[accS2.getSequence().asCodes().at(j2)][accS1.getSequence().asCodes().at(j1)]	// type BP (j2,j1)
							, accS1.getSequence().asCodes().at(i1+1)
							, accS2.getSequence().asCodes().at(i2+1)
							, accS1.getSequence().asCodes().at(j1-1)
							, accS2.getSequence().asCodes().at(j2-1)
							, foldParams))
				;
	} else {
		return E_INF;
	}
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getE_danglingLeft( const size_t i1, const size_t i2 ) const
{
	// Vienna RNA : dangling end contribution
	return Evrna_2_E(vrna_E_ext_stem( BP_pair[accS1.getSequence().asCodes().at(i1)][accS2.getSequence().asCodes().at(i2)]
							  , ( i1==0 ? -1 : accS1.getSequence().asCodes().at(i1-1) )
							  , ( i2==0 ? -1 : accS2.getSequence().asCodes().at(i2-1) )
							  , foldParams
							  ))
			// substract closing penalty
			- getE_endLeft(i1,i2);
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getE_danglingRight( const size_t j1, const size_t j2 ) const
{
	// Vienna RNA : dangling end contribution (reverse base pair to be sequence end conform)
	return Evrna_2_E(vrna_E_ext_stem( BP_pair[accS2.getSequence().asCodes().at(j2)][accS1.getSequence().asCodes().at(j1)]
							  , ( j2+1>=accS2.getSequence().size() ? -1 : accS2.getSequence().asCodes().at(j2+1) )
							  , ( j1+1>=accS1.getSequence().size() ? -1 : accS1.getSequence().asCodes().at(j1+1) )
							  , foldParams
							  ))
			// substract closing penalty
			- getE_endRight(j1,j2);
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getES1( const size_t i1, const size_t j1 ) const
{
#if INTARNA_IN_DEBUG_MODE
	// sanity check
	if (i1>j1) throw std::runtime_error("InteractionEnergy::getES1(i1="+toString(i1)+" > j1="+toString(j1));
	if (j1>=size1()) throw std::runtime_error("InteractionEnergy::getES1() : j1="+toString(j1)+" >= size1()="+toString(size1()));
	if (esValues1 == NULL) throw std::runtime_error("InteractionEnergy::getES1() : ES values not initialized");
#endif

	// return computed value
	return (*esValues1)(i1,j1);
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getES2( const size_t i2, const size_t j2 ) const
{
#if INTARNA_IN_DEBUG_MODE
	// sanity check
	if (i2>j2) throw std::runtime_error("InteractionEnergy::getES2(i2="+toString(i2)+" > j2="+toString(j2));
	if (j2>=size2()) throw std::runtime_error("InteractionEnergy::getES2() : j2="+toString(j2)+" >= size2()="+toString(size2()));
	if (esValues2 == NULL) throw std::runtime_error("InteractionEnergy::getES2() : ES values not initialized");
#endif

	// return computed value
	return (*esValues2)(i2,j2);
}

////////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergyVrna::
isGC( const size_t i1, const size_t i2 ) const
{
	const int bpType = BP_pair[accS1.getSequence().asCodes().at(i1)][accS2.getSequence().asCodes().at(i2)];
	return (bpType==bpCG || bpType==bpGC);
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getE_multiUnpaired( const size_t numUnpaired ) const
{
	return E_type(numUnpaired) * (Evrna_2_E(foldParams->MLbase));
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getE_multiHelix( const size_t j1, const size_t j2 ) const
{
	return (Evrna_2_E(foldParams->MLintern[
	                                     BP_pair[accS2.getSequence().asCodes().at(j2)]
	                                             [accS1.getSequence().asCodes().at(j1)]
	                                    ]));
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getE_multiClosing() const
{
	return (Evrna_2_E(foldParams->MLclosing));
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getEall1() const
{
	// compute Z if needed
	if (E_isINF(Eall1)) {
		Eall1 = computeIntraEall( accS1 );
	}
	return Eall1;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyVrna::
getEall2() const
{
	// compute Z if needed
	if (E_isINF(Eall2)) {
		Eall2 = computeIntraEall( accS2.getAccessibilityOrigin() );
	}
	return Eall2;
}

////////////////////////////////////////////////////////////////////////////


} // namespace

#endif /* INTERACTIONENERGYVIENNA_H_ */
