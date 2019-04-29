# CSparse.NET
A concise library for solving sparse linear systems with direct methods. The code is a C# port of CSparse, written by Tim Davis and part of the [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) project. 

[![Build status](https://ci.appveyor.com/api/projects/status/ji9ig6jrh8u45tb8?svg=true)](https://ci.appveyor.com/project/wo80/csparse-net)
[![Nuget downloads](http://wo80.bplaced.net/php/badges/nuget-dt-csparse.svg)](https://www.nuget.org/packages/CSparse)

## Features

* Sparse LU, Cholesky and QR decomposition of real and complex systems
* Fill-reducing orderings
* Dulmage-Mendelsohn decomposition

All methods are described in detail in the excellent textbook _Direct Methods for Sparse Linear Systems, SIAM, Philadelphia, PA, 2006_ by Tim Davis.

## Examples

* Creating a [sparse LU factorization](https://github.com/wo80/CSparse.NET/wiki/Sparse-LU-example)
* Creating a [sparse Cholesky factorization](https://github.com/wo80/CSparse.NET/wiki/Sparse-Cholesky-example)
* Creating a [sparse QR factorization](https://github.com/wo80/CSparse.NET/wiki/Sparse-QR-example)
* Using [Math.NET Numerics and CSparse.NET](https://github.com/wo80/CSparse.NET/wiki/Math.NET-Numerics-and-CSparse)

## License

    CSparse: a Concise Sparse Matrix package.

    Copyright (c) 2006-2017, Timothy A. Davis.

    This library is free software; you can redistribute it and/or modify it under the
    terms of the GNU Lesser General Public License as published by the Free Software
    the Free Software Foundation; either version 2.1 of the License, or (at your option)
    any later version.

    This library is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along
    with this library; if not, write to the Free Software Foundation, Inc., 51 Franklin
    Street, Fifth Floor, Boston, MA  02110-1301 USA
