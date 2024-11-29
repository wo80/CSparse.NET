### Version 4.2.0 - 2024-09-15

* Make `SymbolicColumnStorage` class public and update `StronglyConnectedComponents` and `DulmageMendelsohn` decomposition accordingly.

### Version 4.1.0 - 2024-06-14

* Add overload for creating a sparse matrix from an enumerable of `ValueTuple`.
* Add matrix `EnumerateIndexedAsValueTuples()` to enumerate entries as `ValueTuple`.

### Version 4.0.0 - 2024-04-03

The major version change is due to the removal of obsolete methods in the `Converter` class. Visibility of that class was changed from public to internal. In case those obsolete methods were still used, please switch to the static conversion methods provided by the `SparseMatrix` class.

Additional changes:

* Add helper method `Helper.ValidateStorage(...)` to validate the structure of a sparse matrix.
* Update to `GetHashCode()` method of `CompressedColumnStorage` class.
* Improvements to documentation.

### Version 3.8.1 - 2023-11-15

* Add overloads for permutation `Invert()` and `IsValid()` methods taking the permutation length as argument.

### Version 3.8.0 - 2023-05-20

* Add overloads for the factorization `Solve()` methods taking `Span<T>` as argument. Note that this introduces a dependency on `System.Memory` for the netstandard2.0 assembly.

### Version 3.7.0 - 2022-05-04

* Add sparse matrix `OfDiagonals` static method (similar to MATLAB spdiags).

### Version 3.6.0 - 2021-11-25

* Remove .NET 4.5 target framework, upgrade .NET 5.0 to 6.0.
* Add constructor that takes explicit non-zeros count to `CoordinateStorage` class.

### Version 3.5.0 - 2021-01-14

* Remove .NET 4.0 target framework, add .NET 5.0.

### Version 3.4.9 - 2020-11-06

* Add `CoordinateStorage` constructor that uses existing storage arrays.
* Convert `CoordinateStorage` to sparse matrix in place.

### Version 3.4.7 - 2020-08-28

* BREAKING: make `SparseLDL` constructor private (use static create methods instead).
* Add complex version of `SparseLDL`.
* Add `matrix.EnumerateIndexed(action)` overload.

### Version 3.4.6 - 2020-07-21

* Add `SolveTranspose` method for `SparseQR`.

### Version 3.4.5 - 2020-06-11

This release introduces the static `SparseMatrix.AutoTrimStorage` option, which enables control over hidden memory allocations in matrix addition and multiplication. By default, the matrix storage will be resized to exactly fit the non-zeros count, which involves new memory allocations. If you want to avoid this, set `AutoTrimStorage` to `false`.

Additional changes:
	
* Add public helper methods `Helper.TrimStrorage(...)` and `Helper.SortIndices(...)`
* Add `DenseMatrix.OfJaggedArray(...)`

### Version 3.4.3 - 2020-05-25

* Add a sparse matrix multiplication overload that accepts the result matrix as a parameter.

### Version 3.4.2 - 2020-05-13

* Make CSparse.NET CLS compliant
* Mark public methods of Converter class as obsolete

### Version 3.4.1 - 2019-10-02

* Improved validation of matrix constructor arguments
* Fixes an issue with `CoordinateStorage` throwing `IndexOutOfRangeException` (introduced in v3.4.0)

### Version 3.4.0 - 2019-09-15

* Parallel dense and sparse matrix multiplication (by Andreas Girgensohn)
* General performance improvements for sparse matrix addition and multiplication

### Version 3.3.0 - 2019-04-29

* Support more target frameworks (including netstandard2.0).
* Public access to members of Dulmage-Mendelsohn decomposition.
* Compute strongly connected components.

### Version 3.2.3 - 2018-11-30

* Added matrix creation helper (e.g. call `SparseMatrix.OfIndexed(s)` to convert coordinate storage).

### Version 3.2.2 - 2018-10-12

* Added MatrixMarket writer.
* BREAKING: make `IProgress interface` compatible with .NET 4.5.

### Version 3.2.1 - 2018-09-17

### Version 3.2.0 - 2018-03-09

* Added new `DenseMatrix` type.
* BREAKING: removed deprecated `CompressedColumnStorage` type.
* BREAKING: removed deprecated `matrix.Norm(int)` method.

### Version 3.1.10 - 2018-03-06

* Rename `CompressedColumnStorage` to `SparseMatrix`.
* BREAKING: `matrix.Multiply(x, y)` overwrites y (instead of update).
* BREAKING: sparse matrix `PermuteColumns` returns a new matrix (instead of update).

### Version 3.1.9 - 2017-01-06

* BREAKING: use static `Create` methods (e.g. `SparseLU.Create(...)`) instead of constructors.

### Version 3.1.4 - 2015-09-19

* Initial release of CSparse.NET (based on Tim Davis CSparse version 3.1.4)
