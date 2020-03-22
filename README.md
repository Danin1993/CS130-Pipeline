# CS130-Pipeline
Graphics Pipeline assignment for Dr. Shinar's CS 130 class (Winter 2020).

## To Get Started
The ideal place to run this project is on any lab machine in Winston Chung Hall at UCR. If you wish to run this at home, you'll need [SCons](https://scons.org/pages/download.html). Clone the project, open up the **CS130-Pipeline** directory, and run `scons`. The project should compile from there.

## Running the Project
This project takes input from 00.txt - 29.txt and renders a scene based on that information. To render one image, run `./driver -i [##].txt -s [##].txt`, where "[##].txt" is any of the input files provided with the project. For example, to render the first image, run `./driver -i 00.txt -s 00.txt`. This will dump the output to **output.png** where you can see the difference between the input and the scene you rendered. It will also output the level of difference between the two images (if any). To run the grading script, run `./grading-script.py .` to get the total differences for all the images. A pass indicates the two images are similar (as determined by staying within an acceptable margin of error). 
