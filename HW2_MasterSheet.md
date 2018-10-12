#HW2 Rhaodes MasterSheet#


##Question 1##

Please generate a new directory named **foo**.
Next move in to the newly created **foo** directory and create a file called **bar**.
Add your username to the **bar** file and print the contents of the file **bar**.
Remove the **foo** directory and its contents.   

**Answer Question 1**

```
mkdir foo
cd foo
touch bar
echo  $USER >> bar
less bar
cd ../
rm -r ./foo
```

##Question 2##

Move to an R environment.

```
R
```

```
mymatrix <- matrix(data = 1:12, nrow = 3, ncol = 4)
mymatrix
```
