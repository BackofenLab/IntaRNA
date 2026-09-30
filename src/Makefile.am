
# sub directories to check for Makefiles
SUBDIRS = . IntaRNA bin

################################################################
# extra files to be added to distribution
################################################################

EXTRA_DIST = easylogging++.LICENCE.txt

###############################################################################
# THE INTARNA LIBRARY
###############################################################################

# the library names to build (note we are building static libs only)
lib_LIBRARIES = libeasylogging.a

# where to install the headers on the system
libeasylogging_adir = $(includedir)

# the list of header files that belong to the library (to be installed later)
libeasylogging_a_HEADERS = easylogging++.h

# the sources to add to the library and to add to the source distribution
libeasylogging_a_SOURCES = easylogging++.cc

###############################################################################


# Bundled Kokkos mdspan headers. Keep the upstream directory layout private to
# IntaRNA's installed headers, with the license beside each header tree.
kokkosdir = $(includedir)/IntaRNA
nobase_kokkos_HEADERS = \
	experimental/__p0009_bits/compressed_pair.hpp \
	experimental/__p0009_bits/config.hpp \
	experimental/__p0009_bits/default_accessor.hpp \
	experimental/__p0009_bits/dynamic_extent.hpp \
	experimental/__p0009_bits/extents.hpp \
	experimental/__p0009_bits/full_extent_t.hpp \
	experimental/__p0009_bits/layout_left.hpp \
	experimental/__p0009_bits/layout_right.hpp \
	experimental/__p0009_bits/layout_stride.hpp \
	experimental/__p0009_bits/macros.hpp \
	experimental/__p0009_bits/mdspan.hpp \
	experimental/__p0009_bits/no_unique_address.hpp \
	experimental/__p0009_bits/trait_backports.hpp \
	experimental/__p0009_bits/type_list.hpp \
	experimental/__p0009_bits/utility.hpp \
	experimental/__p1684_bits/mdarray.hpp \
	experimental/__p2389_bits/dims.hpp \
	experimental/__p2630_bits/constant_wrapper.hpp \
	experimental/__p2630_bits/integral_constant_like.hpp \
	experimental/__p2630_bits/strided_slice.hpp \
	experimental/__p2630_bits/submdspan.hpp \
	experimental/__p2630_bits/submdspan_canonicalize_slices.hpp \
	experimental/__p2630_bits/submdspan_extents.hpp \
	experimental/__p2630_bits/submdspan_mapping.hpp \
	experimental/__p2642_bits/layout_padded.hpp \
	experimental/__p2642_bits/layout_padded_fwd.hpp \
	experimental/mdarray \
	experimental/mdspan \
	mdspan/mdarray.hpp \
	mdspan/mdspan.hpp

nobase_kokkos_DATA = mdspan/LICENSE experimental/LICENSE
EXTRA_DIST += $(nobase_kokkos_DATA)
