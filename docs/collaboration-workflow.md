# Collaborative Workflow

![open_source_collaboration_workflow.png](http://jdimatteo.github.io/images/open_source_collaboration_workflow.png)

### Objectives

* stable and usable master branch
* code review and testing so that master branch code is understandable and of high quality
* someone other than original developer verifying test plan passes
* as simple a a process as could possibly work (with the understanding that it can be made more complicated as needed)

### Example Hello World github usage

Bar will update hello.py, and Foo will code review it.

1. issue created with objective
    * e.g. Foo opens the [HelloWorld github repo page](https://github.com/FooBar11112222/HelloWorld), clicks [Issues](https://github.com/FooBar11112222/HelloWorld/issues), and then [New Issue](https://github.com/FooBar11112222/HelloWorld/issues/new)
![NewIssue.png](http://jdimatteo.github.io/images/NewIssue.png)
2. issue assignee creates feature branch or fork from master (branch preferred)
    * e.g. Bar [creates a branch](https://github.com/blog/1377-create-and-delete-branches)
    * ![CreateBranch.png](http://jdimatteo.github.io/images/CreateBranch.png)
3. objective completed in branch
    * clone the branch, e.g. `$ git clone -b greet-user https://github.com/FooBar11112222/HelloWorld.git`
    * develop code, pushing regularly to branch, e.g. `$ git commit -a -m "Added name prompt, for #1"; git push`
        * note that if you [mention](https://github.com/blog/957-introducing-issue-mentions) the issue number followed by a # in the commit message like above, a comment will automatically be added to the issue with a link to the commit
    * add test plan, verify test plan passing
    * ![TestPlanAdded.png](http://jdimatteo.github.io/images/TestPlanAdded.png)
    * ideally also add automated test so that it will continue to automatically be tested in the future (e.g. a diff test, like [the Bam_Liquidator_Bin_Counter_Test Jenkins test](http://tod.dfci.harvard.edu/jenkins/job/Bam_Liquidator_Bin_Counter_Test/), which runs automatically on every commit, and verifies that the counts have not changed since a prior run's calculated counts for a specific bam file)
4. merge request created
    * click the green "Compare & review" button
    * ![PullRequestButton.png](http://jdimatteo.github.io/images/PullRequestButton.png)
    * click the "Click to create a pull request for this comparison" header button
    * ![CreatePullRequestHeaderButton.png](http://jdimatteo.github.io/images/CreatePullRequestHeaderButton.png)
    * Update the pull request as you see fit, then click the "Send pull request" button
    * ![SendPullRequestButton.png](http://jdimatteo.github.io/images/SendPullRequestButton.png)
5. peer code review
    * Foo opens the pull request and reviews the code changes under the "Files Changed" tab
    * ![FilesChanged.png](http://jdimatteo.github.io/images/FilesChanged.png)
    * Foo adds a comment to the diff:
    * ![Feedback.png](http://jdimatteo.github.io/images/Feedback.png)
6. assignee responds to peer comments
    * Bar reads the reviewer comments, and pushes commits to the greet-user branch
    * github automatically adds those commits to the pull request
    * Foo reviews the changes, and signs off on it
    * ![GoAhead.png](http://jdimatteo.github.io/images/GoAhead.png)
7. assignee merges branch with master
    * Bar clicks the "Merge pull request" button
    * ![MergePullRequest.png](http://jdimatteo.github.io/images/MergePullRequest.png)
        * if the pull request can't be merged automatically due to conflicts, read [documentation](https://help.github.com/articles/resolving-a-merge-conflict-from-the-command-line/) or ask another developer for help
    * Bar clicks the "Confirm merge" button
    * the pull request is now marked as merged
    * Bar clicks the "Delete branch" button
        * Note that all branch history is preserved, deleting it just indicates that it is no longer being developed
    * Bar can go ahead and close the issue
8. Jenkins GUI wrapper jobs (and anything else using master) automatically uses updated master

### Prior Drafts

* initial diagram, more complicated than current draft, with stable master but named release branch for a dedicated testing phase:
![open_source_collaboration_workflow_with_release_phase.png](http://jdimatteo.github.io/images/open_source_collaboration_workflow_with_release_phase.png)

### References

* code reviews: http://blog.codinghorror.com/code-reviews-just-do-it/
* git:
    * http://scottchacon.com/2011/08/31/github-flow.html
    * https://github.com/blog/1557-github-flow-in-the-browser
    * https://help.github.com/articles/using-pull-requests
    * https://github.com/blog/1377-create-and-delete-branches
    * https://github.com/blog/957-introducing-issue-mentions
    * https://help.github.com/articles/how-do-i-set-up-a-team
    * http://nvie.com/posts/a-successful-git-branching-model/
    * https://help.github.com/articles/resolving-a-merge-conflict-from-the-command-line/
    * https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging#Basic-Merge-Conflicts
    * [Version Control with Git](http://www.amazon.com/Version-Control-Git-collaborative-development/dp/1449316387)
