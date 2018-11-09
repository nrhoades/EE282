<h1 id="toc_0">HW2 Rhoades MasterSheet</h1>

<h2 id="toc_1">Question 1</h2>

<p>Please generate a new directory named <strong>foo</strong>.
Next move in to the newly created <strong>foo</strong> directory and create a file called <strong>bar</strong>.
Add your username to the <strong>bar</strong> file and print the contents of the file <strong>bar</strong>.
Remove the <strong>foo</strong> directory and its contents.   </p>

<p><strong>Answer Question 1</strong></p>

<div><pre><code class="language-none">mkdir foo
cd foo
touch bar
echo  $USER &gt;&gt; bar
less bar
cd ../
rm -r ./foo</code></pre></div>

### Question 1 Comments:
Very well done thank you. Also just an fyi, if you want to avoid getting prompts when deleting you can use the switch -rf intead of -r, just be very careful when you do this and don't delete all your data, especially when using wildcards i.e rm -rf * <-- very dangerous.

<h2 id="toc_2">Question 2</h2>

<p>Move to an R environment.
We will be using the built in Dataframe <em>mtcars</em> for this assignment.</p>

<div><pre><code class="language-none">R</code></pre></div>

<p>What is the name of the index of the first column of <em>mtcars</em>.</p>

<p>Which of the following syntax produce the same output?
Are any of the syntax unable to access the information in the dataframe?
(mtcars[,&#39;mpg&#39;] vs. mtcars[&#39;mpg&#39;] vs. mtcars$mpg vs. mtcars[[&#39;mpg&#39;]])</p>

<div><pre><code class="language-none">mtcars[,&#39;mpg&#39;]
mtcars[&#39;mpg&#39;]
mtcars$mpg
mtcars[[&#39;mpg&#39;]]</code></pre></div>

<p><em>ANSWER</em><br>
mtcars$mpg<br>
mtcars[[&#39;mpg&#39;]<br>
mtcars[,&#39;mpg&#39;]<br>
All return a list of values from the from the MPG column without the corresponding keys.</p>

<p>While mtcars[&#39;mpg&#39;] returns the list of values with the cars that correspond to them.</p>

<hr>

<p>Convert the the dataframe to a matrix and reassess tun mtMatrix[,&#39;mpg&#39;].</p>

<div><pre><code class="language-none">mtmatrix &lt;- as.matrix(mtcars)
mtmatrix [,&#39;mpg&#39;]</code></pre></div>

<p>How do your results change?</p>

<p><em>ANSWER</em><br>
You are still able to access the same data in the matrix, but each mpg is pecfically connected to a car.
Prints differently.</p>

<h2 id="toc_3">Question 3</h2>

<p>Outside of an R environmnt</p>

<p>Create 3 file titled <em>foo.sh</em></p>

<p>Write a script to print &quot;Hello World&quot; in the terminal.</p>

<p>Using Octal permissions, make files respectively executable by everyone, but only writeable by the owner.</p>

<p>Execute the <em>foo.sh</em> script.</p>

<h4 id="toc_4">Answer Question 3</h4>

<div><pre><code class="language-none">touch foo.sh
echo  -e &#39;#!/bin/bash/\necho Hello World&#39;  &gt;&gt; foo.sh
chmod 744 foo.sh
sh foo.sh</code></pre></div>
