# List of things to check before submitting the code for review

This document presents the common mistakes often seen during the review and
also includes some general guidelines in code writing and testing.

Please make sure as a PR (Pull-Request) submitter
(or as a first-stage code reviewer) these are eliminated wherever possible.

## General

-   Use PR / commit descriptions where possible to indicate what has
    changed, it helps your peer reviewer to understand basic differences
    among 20+ commits

-   Don't change the status of the associated Jira ticket(s to implemented unless merged and
    **agreed with your reviewer**

-   Include ticket number in the Pull-Request message (the PR template will
    have an entry for you to fill in)

## Builds

-   Don't submit a branch for a review unless it passes all pre-submit CI
    jobs

-   Check that your code didn't add any new compiler warnings

## Library code

-   Don't reinvent the *wheel*, reuse existing infrastructure

    Chance that you need something slightly different which tempts
    you to write a very similar function to an existing one is
    likely an indicator that you should think twice.

-   **Do NOT copy/paste code**. There might be rare occasions where it is
    useful but in general it denotes bad code design (particularly for
    bigger code blocks). Consider using a function or a template
    function with parameters.

-   Calls to functions (e.g., internal ones) should check for their
    return `status` which is propagated down the call stack in transparent
    manner (use as appropiate `return error(...);` or `return status;`)

-   Use variable names matching their common interpretation and ideally
    only in that context

    e.g., `i` ~ row index, `j` ~ column index

-   Use `const` keyword on pointers and references wherever required
    (including the internal functions) but don't use them for
    scalars passed-by-value.

    Use `const` also when casting constant pointers (otherwise you might wipe
    this protection off).

-   Document your code, particularly pay extra attention to

    -   New data structures, enums, etc.

    -   Internal functions with many parameters: is it obvious how all
        parameters are connected? Will it ever get reused?

-   Test all aspects of your input early in your user-facing API

    -   Never dereference a pointer without checking earlier that it is
        not `nullptr` (or `NULL`)

## Documentation

-   Build the documentation and read it from there if the formatting is
    correct

    -   Is the API listed in `.rst` files that the documentation appears
        at all?

-   Imagine that you read the documentation without prior knowledge of
    the API, would you be able to use it? Write examples? Understand
    what needs to be on input and how to process output?

-   Are all errors the API can return listed?

-   Follow the common patterns how to refer to things, examples:

    -   If your API has an example, include it in the documentation via *collapsed
        file loading*

-   Typography: don't place **spaces** before commas, semicolons,
    colons etc. (e.g., not *~~word1 , word2 : word3~~* but *word1,
    word2: word3*)

## Example

-   Check returned status of all public APIs called

-   Include a basic check if the API works and return non-zero code if
    it didn't, this helps to use the example as a test in the CTest framework

## Testing

-   Good practice is to run only your new unit tests and check the
    coverage. Do you see what you expected? Ideally coverage LOC percentage would
    not decrease

-   When creating wrong input cases, make sure that exactly one thing is
    wrong (e.g., if matrix is supposed to be symmetric and of matching
    dimensions, there should be multiple test cases such as the
    dimension is right but it is not symmetric, it is symmetric but the
    dimension is wrong)

    -   Also make it just slightly wrong, not by miles (e.g., if we
        expect $n=m$, choosing 5 & 6 might be more likely to fail than 5
        & 234891)

-   Choose such input that it is likely to detect error, examples what
    not to do:

    -   square matrix $m=n$ so making a mistake by swapping $n$ & $m$ is not
        discovered

    -   input vector or matrix of all ones doesn't test that the correct indices are
        accessed

    -   Complex numbers with same real and imaginary parts

-   Test on edge case input as well as typical

    -   e.g., 0-sizes, 0-nnz, 1 row matrix, ...

-   Set correct tolerances for the expected results, use tolerance relative to machine epsilon

-   Use standardized facilities to check if your results match with the
    reference, **don't write comparison functions that already exist**

-   If using radnomly generated data, initialize the random seed to allow for
    reproducibility

-   In positive quick exit tests, check the validity of the output
    (e.g., that the matrix generated which is supposed to be empty is
    indeed valid, of the right dimension and empty)

-   It might not be necessary to test a complete combination of all
    parameters with the same input, although sometimes it might be desirable. It
    might make more sense to vary the input for different types to get
    better coverage with smaller number of tests.
