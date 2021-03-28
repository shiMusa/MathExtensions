# MathExtensions

## Advanced Math functions for Jai

The goal of this project is to progressively write a comprehensive
mathematical library for Jai, similar to the Gnu Scientific Library or Java Apache Math.

The focus lies in providing the functionality first, performance
second. 

Feel free to add functionality and performance upgrades.

## Warning

I'm not a very good programmer. I make mistakes, write unoptimal code, make wrong design decisions.
As a scientist, I'm an end-user of these math libraries. My goal is to make this library simple to use and easy to understand.

## Vision

Since Jai is as low-level as it gets, the utopia would be to write ALL of the math in Jai.

However, there are decades worth of work already put into libraries such as BLAS, LAPACK, etc. So maybe wrappers for those libraries are easier and more stable and reliable than writing everything from scratch.

I'll start by writing naive Jai implementations and later down the line we can improve incrementally.


## Structure

Let's keep it as simple as possible:

- Operator overloading where it makes sense (linear algebra, complex numbers, quaternions, etc.), simple function calls otherwise.
- No special type-aliases if not stricly necessary
- A file for each major topic.

## Important

Give sources for algorithms written so that we can take a look and help debugging.

Write (at least some) tests contained in each file.

## TODOs

A list of open tasks for the categories already working on:

- Flags for matrices (anti-/symmetric, hermitian, upper/lower triangle, diagonal, unitary etc.) for using specialized algorithms improving performance
