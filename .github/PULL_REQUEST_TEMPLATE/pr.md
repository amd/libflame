Peer Review Checklist for Pull-Request (PRPR)
=============================================

*Revision 19/6/2025*

See Appendix for the *Peer Review Guide*.

**Synopsis** <fill in>

JIRA Tickets #<fill in>

Fixes #<fill in>

Relates to #<fill in>

Reviewers: Add reviewers in the right-hand-side tab, select one.

# Self-Review Checklist

This checklist is used to help the developer(s) ensure that the new code
is ready to be peer reviewed and the Pull-Request issued.

**Before passing to the PRs, this checklist needs to be completed**. If any
answers are not `Yes` or `NA`, please provide details to the reviewer so that any
problems get resolved.

**Don't waste PR's time by not checking this basic Self-Review**

<table border="1">
    <tr>
        <th>YES</th>
        <th>NA</th>
        <th>Item</th>
        <th>Description</th>
    </tr>


<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>1</td> <td> Does the Pull-Request message reference the correct JIRA ticket(s)?  </td> </tr>

<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>2</td> <td>Does the commit message reference PR issue number?  </td> </tr>

<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>3</td>
<td>

Is the PR from a branch branched off `amd-main` and follows the required convention  `u/<ntid>/PR/<branch name>`?

 * Are there any conflicts to be resolved? Merge/Rebase and resolve as required.
 * Are all the automatic builds of the branch successful? Fix/update as required.

</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>4</td>
<td>

Is the code polished with `clang-format`?<br>

</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>5</td>
<td>
Is the standard copyright banner and license information at the top of all relevant files?<br>
Also, are/is the copyright year(s) on the banners correct? Copyright is expressed on a
continuous range of years since its inception, so must have at least (e.g. 2023)
and optionally followed by the last year the file was updated (e.g. 2023-2025).
</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>6</td>
<td>

Does the **full-debug build** complete successfully without any "legit" warnings?<br>
**Add only "legit" warnings to the Suppressions List.**

</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>7</td>
<td>
Has the C and C++ API been sufficiently unit-tested?<br>

 * Code coverage (lines of code)
 * Algorithmic paths coverage (refinement of a solution under bad problem conditioning)
 * Ideally at least one unit-test for each error returned by the API
 * Tests for both row-major and column-major data if appropriate

</td>
</tr>
</table>

<br>
Only for new or updated API requiring documentation

- [ ] Ignore next table, section not relevant

<table border="1">
<tr>
<th>YES</th>
<th>NA</th>
<th>Item</th>
<th>Description</th>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>1</td>
<td>

Does the documentation correctly describe the API?
 * A high-level description of the implemented algorithm is provided at the top of the API (i.e., comments used for Doxygen) or where applicable; if applicable, include the maximum expected error or tolerance
 * All inputs and output are sufficiently explained and marked as `in`/`out`
 * List and explain all returned errors (if advice how to rectify the error is available, include it in the description)
 * If there are example sources, refer/link/include these in the documentation

</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>2</td>
<td>Has the API documentation been checked for spelling and grammar?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>3</td>
<td>
Are all references to groups of implementations for the same API using <code>\p api_name_?</code>?<br>
</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>4</td>
<td>Does compiling the documentation, i.e., running Sphinx/Doxygen/LaTeX report no errors?</td>
</tr>
</table>


If all questions are answered YES, then assign reviewer to PR. It is strongly
advised for the developer(s) to also review the PR checklist prior to
contacting the reviewer.

**Proceed with the PR only if the above table(s) is/are completed**

# Peer review (PR) checklist

## Purpose and Scope

This checklist is to assist the peer reviewer in assessing the quality
of the new code (project files) before declaring it to be 'complete'.

The code should not be considered complete **until all items** are answered
`Yes`/`NA` or justified explain the discrepancy. The
PR Checklist is composed of four major parts:

 1. Documentation,
 2. Code,
 3. Examples, and
 4. Unit testing.

All are essential to the review process.

## Documentation (compulsory for new and updated APIs)

- [ ] Ignore next table, section not relevant

<table border="1">
    <tr>
        <th>YES</th>
        <th>NA</th>
        <th>Item</th>
        <th>Description</th>
    </tr>
    <tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>1</td>
<td>
    Does the documentation correctly describe the API?

* A high-level description of the implemented algorithm is provided at the top of the API (i.e., comments used for Doxygen) or where applicable; if applicable, include the maximum expected error or tolerance
* All inputs and output are sufficiently explained and marked as in/out
* List and explain all returned errors (if advice how to rectify the error is available, include it in the description)
* If useful, refer to the examples
</td>
</tr>
<tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>2</td>
<td>Is the documentation clear, concise, and free of any grammatical or spelling errors?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>3</td>
<td>Does the documentation state clearly where non-successful exit statuses still convey useful or relevant results? Some APIs go beyond success/failure, e.g., in numerical optimization, a warning can be issued to notify a sub-optimal solution.</td>
</tr>
<tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>4</td>
<td>Optional, are any suitable references to further resources provided?</td>
</tr>
</table>

## Code

- [ ] Ignore next table, section not relevant


<table border="1">
<tr>
<th>YES</th>
<th>NA</th>
<th>Item</th>
<th>Description</th>
</tr>
<tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>1</td>
<td>Have the API goals been achieved? (e.g., target performance, etc.)</td>
</tr>
<tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>2</td>
<td>Does the API correctly implement the algorithm(s) (e.g., for a matrix-vector product API, does the code correctly perform a matrix-vector product for all valid inputs)? Does it handle row-major and column-major data?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>3</td>
<td>Is the code reasonably easy to read and understand? Is it clear in the code how the underlying algorithm is implemented? Do you consider that other colleagues will be able to debug the API?

Has the source code been written in an appropriately concise and efficient way? Do you understand what is happening? Speak to the developer if you have any doubts and suggest improvements.</td>
</tr>
<tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>4</td>
<td>Is the source code appropriately documented? Are the comments relevant and helpful? Are descriptions and comments correctly styled and relevant? Do you understand all the comments in the code? If not, then a user is also unlikely to understand them so query any such comments with the author. While comments in the code don't have to be grammatically correct (e.g., single line comments don't require full sentences), they should still make sense.</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>5</td>
<td>Have you been (reasonably) unsuccessful in breaking the code? The main objective here is to assess the robustness of the API. For example, by testing how well it handles unintended usage either by passing unexpected (valid input) data, or by defining unexpected combination of algorithmic options, etc. This requires skillful creativity from the peer-reviewer.</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>6</td>
<td>Does the code gracefully handle trivial and edge-cases (e.g., factorizing identity matrix, working with matrix of dimension 1x1 or sparse matrix with no nonzeroes)?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>7</td>
<td>Is there adequate error handling?</td>
</tr>

<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>8</td>
<td>Is the code thread-safe?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>9</td>
<td>Are variable and function names short, sensible, and consistent with those of other similar functions?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>10</td>
<td>Does the API signature have a <code>const</code> qualifier for all input parameters passed by reference?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>11</td>
<td>Do all the source files have copyright banner with correct year(s) as well as license information?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>12</td>
<td>Is the source code formatted (polished) appropriately?<br>

Also, if there are changes to `CMakeList.txt` files, are these polished?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>13</td>
<td>

Are all code files correctly placed (i.e., in the
`src/` or `include/` directories)?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>14</td>
<td>Have all necessary additional APIs been written? e.g. C interfaces </td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>15</td>
<td>The code does not produce any warnings when compiling with full-debug build using GCC and AOCC (including Address SANitizer checking)?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>16</td>
<td>Do all periodic Jenkins jobs (Linux/Windows) pass without producing significant compiler warnings?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>17</td>
<td>Optional, look closely at the conformance tests. Are there any potential tests that have not been included? This is worth spending some time thinking about.</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>18</td>
<td>Optional, is the new or updated API performance adequately efficient?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>19</td>
<td>Optional, is the parallel implementation (OpenMP) showing acceleration?</td>
</tr>
<tr>
<td>

- [ ]  &nbsp;

</td>
<td>

- [ ]  &nbsp;

</td>
<td>20</td>
<td>

Are all array-type (vectors, etc) that are dynamically allocated ... properly de-allocated 
</td>
</tr>
</table>

## Examples


- [ ] Ignore next table, section not relevant

<table border="1">
    <tr>
        <th>YES</th>
        <th>NA</th>
        <th>Item</th>
        <th>Description</th>
</tr>
<tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>1</td>
<td>Is the example easy to understand, solving a simple (yet not completely trivial) case and can it be used as a template by the user?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>2</td>
<td>Does the example avoid numerical issues (e.g., close to zero division, poor matrix conditioning, etc.)?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;
</td>
<td>

- [ ] &nbsp;
</td>
<td>3</td>
<td>Given the example input data, is the output uniquely defined?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>4</td>
<td>Does the example check the returned status of every relevant API call and gracefully handle the condition?
For example, on non-successful status return, identify if the API still returns relevant results (e.g., iterative methods stop with sub-optimal solutions).</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>5</td>
<td>Does the example return exit status 0 on successful (expected) completion of a run and a non-zero exit status otherwise (CTest compatibility, this allows the use of the examples for testing)?

Also, when printing a solution, print also a message that the result is within the tolerance, where appropriate.</td>
</tr>
</table>

## Unit testing


- [ ] Ignore next table, section not relevant

<table border="1">
    <tr>
        <th>YES</th>
        <th>NA</th>
        <th>Item</th>
        <th>Description</th>
    </tr>
    <tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>1</td>
<td>

Are there tests checking that invalid input is correctly handled by the API returning expected error status? (e.g., negative dimension, `nullptr` or `NULL`  pointer)</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>2</td>
<td>Are edge cases tested (e.g., dimension is 1)?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>3</td>
<td>Are trivial cases tested (e.g., when providing the solution as input)?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>4</td>
<td>

Are all *reasonable* algorithmic failures tested where possible? (e.g., factorizing singular matrix)</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>

<td>5</td>
<td>Do the unit tests cover a range of realistic problems?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>6</td>
<td>Is the code coverage by the unit tests sufficient? Particularly, are all possible code paths and functionality use cases exercised? (e.g., different AVX extensions, types of matrices, row-major and column-major data)</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>7</td>
<td>Is it easy to understand what is being tested from the source and comments?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>8</td>
<td>Do the automatic Linux/Windows builds build successfully without failing? Are all unit tests passing?</td>
</tr>
<tr>
<td>

- [ ] &nbsp;

</td>
<td>

- [ ] &nbsp;

</td>
<td>9</td>
<td>

Optional, are precision-related tolerances defined relative to the type (`float`/`double`) using scaling on `epsilon` (machine precision, taken from the appropriate header file) and not hard-coded constants?
</td>
</tr>
</table>

# Appendix - Quick guide

<details>
    <summary>Click to expand</summary>#


## Motivation, purpose, and scope


The peer review process needs to ensure that all AMD libraries are
developed to the highest standards of quality in order to deliver the
best experience we can for the user. This means that a diverse set of
factors needs to be considered as part of the review, including the
following.

1.  **Code implementation**: design, correctness, efficiency,
    robustness, and usability

2.  Monitoring of **quality**: checks via thorough testing (**unit
    tests**, coverage, etc.)

3.  User experience -- product learning curve

    a.  **Documentation**: API specs, underlying algorithm description,
        theoretical math background, etc.

    b.  **Examples** as templates for users to start with

The development process involves two parties who are jointly responsible
for a successful completion:

**Developer(s)**: assigned the tasks of developing and documenting the
software. Either a single person or a small team.

**Peer Reviewer (PR)**: a **single** designated person helping the
developer(s) make the right decisions along the way ("development
buddy") and ultimately responsible for the final approval confirming the
overall quality of the code and documentation, as well as for the
thorough testing of the API and example(s).

Peer reviewing is a process that starts at the same time as the planning
of the new API or updating an existing one and culminates with the
merging of the feature branch. It is strongly encouraged that
development be done in an open manner and to have periodic review
meetings with the peer reviewer. The process is presented below as a
flow chart and assumes that development is done using a feature branch
based off a trunk and that revisions are done using Gerrit. These are
not mandatory, and the PR process can be adapted to similar setups
(e.g., GitHub). To aid the peer review process two checklists are
presented. These checklists should be considered as templates which can
be customized to suit the requirements and processes of specific
libraries.

## Notes

**Full-debug build** refers to either a GCC or AOCC (Clang) target that
uses -Wall and -Wextra or similar flags and is the default defined by
CMake build type "Debug." This may also include checks with an address
sanitizer facility.

**Suppressions List** is an optional list of source filenames and
associated compiler/linker warnings that can be safely ignored. Warnings
in this list are suppressed when assessing QA.

**All Jenkins jobs** refers to a collection of (automatic) jobs that are
setup to verify code quality using debug builds on both Linux and
Windows platforms. These jobs not only build but also run all the
unit-tests and examples and report compiler, unit-tests and examples
result logs. It possibly also includes jobs that build a code coverage
report and compile the documentation.

**Unit-tests** refers to all test source codes related to the project.
These may be a suite in Google Tests/CTests.
</details>
