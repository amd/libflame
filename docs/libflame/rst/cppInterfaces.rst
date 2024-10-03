..  Copyright (C) 2024, Advanced Micro Devices. All rights reserved.

..  Redistribution and use in source and binary forms, with or without
..  modification, are permitted provided that the following conditions are met:

..  1. Redistributions of source code must retain the above copyright notice,
..  this list of conditions and the following disclaimer.
..  2. Redistributions in binary form must reproduce the above copyright notice,
..  this list of conditions and the following disclaimer in the documentation
..  and/or other materials provided with the distribution.
..  3. Neither the name of the copyright holder nor the names of its
..  contributors may be used to endorse or promote products derived from this
..  software without specific prior written permission.

..  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
..  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
..  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
..  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
..  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
..  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
..  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
..  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
..  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
..  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
..  POSSIBILITY OF SUCH DAMAGE.

.. _cppInterfaces:


CPP Interfaces
-------------------

AOCL-LAPACK provides standard C++ Interface to LAPACK APIs.

* :ref:`LinearSolve`
   * :ref:`cholesky`
   * :ref:`LU_General_matrix_driver`
   * :ref:`LU_Computational`
   * :ref:`LDL`
   * :ref:`LDL_computation`
   * :ref:`Triangular_Computational`
   * :ref:`Auxilliary_routines_LinearSolve`
* :ref:`LeastSquare`
   * :ref:`Standard_Least_square`
   * :ref:`Constrained_least_squares`
   * :ref:`Auxiliary_routines_least_square`
* :ref:`SingularValueDecomposition`
   * :ref:`SVD_Computational_Routines`
* :ref:`EigenvalueRoutines`
   * :ref:`Non_symmetric_eigenvalues`
   * :ref:`Symmetric_eigenvalues`
* :ref:`OrthogonalUnitaryFactors`
   * :ref:`QR`
   * :ref:`LQ`
   * :ref:`CosineSineDecomposition`
   * :ref:`HouseholdReflector`
   * :ref:`GivensJacobiPlaneRotations`
* :ref:`AuxilliaryRoutines`
   * :ref:`Parameters`
* :ref:`BLASLike`
   * :ref:`InitializeCopyConvert`
   * :ref:`MatrixNorm`
   * :ref:`ScalarOperations`
   * :ref:`Level1BLASLike`
   * :ref:`Level2BLASLike`
   * :ref:`Level3BLASLike`


.. toctree::
   :hidden:
   :caption: Contents:

   LinearSolve.rst
   LeastSquare.rst
   SingularValueDecomposition.rst
   EigenvalueRoutines.rst
   OrthogonalUnitaryFactors.rst
   AuxilliaryRoutines.rst
   BLASLike.rst