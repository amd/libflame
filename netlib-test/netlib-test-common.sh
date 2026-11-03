#!/bin/bash

#
# Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
#

# Select an unversioned shared, versioned shared, or static library in that order.
select_default_library()
{
	local directory=$1
	local base_name=$2
	local static_name="${base_name}.a"
	local shared_candidates

	# Prefer the conventional unversioned shared library.
	if [[ -f "$directory/${base_name}.so" ]]
	then
		basename "$directory/${base_name}.so"
		return
	fi

	# Accept a versioned shared library when the unversioned name is absent.
	shared_candidates=("$directory/${base_name}.so."*)
	if [[ -f ${shared_candidates[0]} ]]
	then
		basename "${shared_candidates[0]}"
		return
	fi

	# Fall back to the static archive.
	if [[ -f "$directory/$static_name" ]]
	then
		echo "$static_name"
		return
	fi

	return 1
}

# Return an absolute library directory, interpreting relative paths from the copy directory.
resolve_library_directory()
{
	local directory=$1
	local copy_directory=$2

	if [[ $directory = /* ]]
	then
		realpath -m -- "$directory"
	else
		realpath -m -- "$copy_directory/$directory"
	fi
}

# Resolve and validate test libraries, then configure their runtime search path.
configure_test_libraries()
{
	local copy_directory
	copy_directory="$(pwd)/libflame_netlib"

	# Keep directories absolute so later validation, copying, and loader lookup share one path.
	BLAS_LIB_PATH=$(resolve_library_directory "$BLAS_LIB_PATH" "$copy_directory") || return 1
	LAPACK_LIB_PATH=$(resolve_library_directory "$LAPACK_LIB_PATH" "$copy_directory") || return 1

	# Auto-select BLAS when its filename was not provided.
	if [[ -z $BLAS_LIB ]]
	then
		BLAS_LIB=$(select_default_library "$BLAS_LIB_PATH" "libblis-mt") || {
			echo "No libblis-mt shared or static library found in $BLAS_LIB_PATH"
			return 1
		}
	elif [[ ! -f "$BLAS_LIB_PATH/$BLAS_LIB" ]]
	then
		echo "BLAS library $BLAS_LIB_PATH/$BLAS_LIB does not exist"
		return 1
	fi

	# Auto-select LAPACK when its filename was not provided.
	if [[ -z $LAPACK_LIB ]]
	then
		LAPACK_LIB=$(select_default_library "$LAPACK_LIB_PATH" "libflame") || {
			echo "No libflame shared or static library found in $LAPACK_LIB_PATH"
			return 1
		}
	elif [[ ! -f "$LAPACK_LIB_PATH/$LAPACK_LIB" ]]
	then
		echo "LAPACK library $LAPACK_LIB_PATH/$LAPACK_LIB does not exist"
		return 1
	fi

	# Report the resolved test configuration.
	echo "BLAS_LIB_PATH = $BLAS_LIB_PATH"
	echo "BLAS_LIB = $BLAS_LIB"
	echo "LAPACK_LIB_PATH = $LAPACK_LIB_PATH"
	echo "LAPACK_LIB = $LAPACK_LIB"
	echo "LAPACK_TEST_DIR = $LAPACK_TEST_DIR"
	echo
	echo
	echo "**********************************"
	echo

	# Add library directories to the runtime loader's search path.
	RUNTIME_LIBRARY_PATH="$LAPACK_LIB_PATH:$BLAS_LIB_PATH"
	export LD_LIBRARY_PATH="$RUNTIME_LIBRARY_PATH${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
}
