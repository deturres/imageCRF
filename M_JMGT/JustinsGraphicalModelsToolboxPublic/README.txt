DOCUMENTATION

http://phd.gccis.rit.edu/justindomke/JGMT/


QUICK START

1) Go to the root directory in matlab

2) Compile the code.  You shouldn't need to do this if binaries are available for your system (64 bit linux or 64 bit os X).  If you do need to compile the binaries, please email them the author so they can be included for others.

>> mex -setup % make sure you have a c++ compiler installed
>> compile

3) (optional) Compile the multi-threaded code back-TRW code.  You probably _will_ need to do this if you want to use the multi-threaded code, because binaries using openmp don't appear to be very portable.  Your compiler must support openmp in order to do this.  Recent versions of gcc, icc, visual studio all do, though this has only been tested using gcc

>> compile_openmp

4) Add the libraries to the path

>> addpath(genpath('.'))

5) Run an example of training and evaluating a CRF

>> example_binarydenoising

6) A more complex example is available from Examples/example_backgrounds  (requires downloading some code)


DIFFERENT WAYS YOU COULD USE THE TOOLBOX

1 - Just call the [Inference] methods (trw/meanfield), do everything else on your own.

2 - Use the [Differentiation] methods (bprob_trw, implicit_diff) to calculate parameter gradients by providing your own loss functions.  Do everything else on your own.

3 - Use the [Loss] methods (em, implicit_loss) to calculate parameter gradients by providing a true vector x and a loss name

4 - Use the [CRF] methods to calculate calculate almost everything (deal with parameter ties for a specific type of model, etc.)


QUESTIONS

justin.domke@rit.edu