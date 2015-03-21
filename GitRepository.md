The branches are divided with respect to the level of stability of the igatools capabilities.
  * **Experimental**: code that is currently under development and possibly without documentation. _Interfaces will change almost surely_. Branches in the **Experimental** state are candidates to be integrated into one (or more) **Unstable** branch. The branches currently in this state are:
    * **hierarchical**: branch to be used for the development of hierarchical spaces
    * **multipatch**: branch to be used for the development of the multipatch support of igatools.
    * **parallel**: branch to be used for the development of the parallel support of igatools (Intel TBB and MPI).
    * **linear\_algebra**: branch to be used for the development of linear-algebra support of igatools and wrappers to external linear algebra packages (e.g. Trilinos, PETSc).

  * **Unstable**: code that is tested and partially documented. Interfaces may change in the future. Branches in the **Unstable** state are candidates to be integrated into the **Stable** branch. The branch currently in this state is:
    * **develop**: this branch contains the version of igatools that is under the supervision of the leading developers.

  * **Stable**: code that is tested, cleaned and documented. Interfaces will not change (hopefully). No branches are present at this time.


We encourage to work with local branches as much as possible and then push your commits to one of the Experimental branches. We ask to follow some (simple) rules:
  * in Google Code we cannot assign different permissions to different branches, so **_DON'T PUSH YOUR COMMITS TO UNSTABLE OR STABLE BRANCHES!!!_**. Please ask the **leading developers** to merge one or more Experimental branch into an Unstable branch. Any violation of this rule will be discussed by the **leading developers** that can decide to revoke your permission to push commits to the igatools project.
  * pushed commits to an Experimental branch should not break the library, i.e. all pushed commits should work (or at least should let the library to be compiled and used without the new feature). The reason for this is the fact that other people may work on the same Experimental branch. If you need to create sub-branches, do it locally on your machine or ask a leading developer to do it for you in Google Code.
  * We advocate the [Test Driven Development model](http://en.wikipedia.org/wiki/Test-driven_development), so any new feature introduced in a Experimental branch should be validated following the [igatools testing policy](Tests.md).
