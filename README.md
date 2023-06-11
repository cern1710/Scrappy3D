# Scrappy3D

Scrappy3D is a 3D modeling, rendering, and animation package based on [15-462/662 Computer Graphics](http://15462.courses.cs.cmu.edu) at Carnegie Mellon University.

## Disclaimer

This is merely for educational purposes and is not intended to be used as reference. The repo has been intentionally renamed to make it harder to find.

## General Setup

1. Install a C++ compiler:
  - Windows: Visual Studio 2022
  - MacOS: XCode (latest available)
  - Linux: g++ (latest available)
2. Install [node](https://nodejs.org) (our build system is written in command-line javascript.)
3. Clone this repository.
4. Download and extract the nest-libs as a child of the repository folder:
  - Linux: https://github.com/15-466/nest-libs/releases/download/v0.10/nest-libs-linux-v0.10.tar.gz
  - MacOS: https://github.com/15-466/nest-libs/releases/download/v0.10/nest-libs-macos-v0.10.tar.gz
  - Windows: https://github.com/15-466/nest-libs/releases/download/v0.10/nest-libs-windows-v0.10.zip


## Building and Running

Run, from a command prompt that has your compiler available:
```
#change to the directory with the repository:
$ cd Scrappy3D
#use Maek to build the code:
Scrappy3D$ node Maekfile.js
#run the UI:
Scrappy3D$ ./Scrappy3D
#run the tests:
Scrappy3D$ ./Scrappy3D --run-tests
```
Note that you _should_ read `Maekfile.js`. about the available command line options and how to configure your own build. All the code has been nicely documented to help you understand the building process and reinforce your learning.

Below is a non-exhaustive list of common build issues along with their suggested solutions. Let us know if you encounter a problem that is not addressed here.

### macOS
> `unknown warning option '-Wno-unused-but-set-variable'`

`'-Wno-unused-but-set-variable'` is a [new compiler flag introduced in Clang 13.0.0](https://releases.llvm.org/13.0.0/tools/clang/docs/ReleaseNotes.html#new-compiler-flags), so you may encounter this error if your Clang is not up to date. To resolve this, do the following:

1. Check that your macOS version is Monterey or above, then upgrade Xcode to the latest version (should be at least Version 13.0).
2. Open Xcode and navigate to the top menu bar. Select **Xcode->Preferences->Locations**, then set **Command Line Tools** to your latest Xcode version (e.g. Xcode 14.2 (14C18)). 
3. Run `clang -v` in your terminal to check if thel Clang version is $\ge$ 13:
```
Apple clang version 14.0.0 (clang-1400.0.29.202)
Target: ...
Thread model: ...
InstalledDir: ...
```
4. Re-try building Scrappy3D. If the error still shows up, go to `Maekfile.js` and comment out the line where the flag is located. We included this flag because it is beneficial to your learning and omitting it may cause it to not build on some other computers.

### Linux
> `/usr/bin/ld: cannot find -lasound`

Your machine is missing a package. Simply install it by running `apt-get install libasound2-dev`.

### Windows
> `library machine type 'x64' conflicts with target machine type 'x86'`

You are probably building from a shell with the developer tools set to target x86 (older, 32-bit instruction set) instead of x64 (modern, 64-bit flavor of x86, also sometimes known as “x86_64”). Check [this guide](https://learn.microsoft.com/en-us/cpp/build/how-to-enable-a-64-bit-visual-cpp-toolset-on-the-command-line?view=msvc-170) for enabling an x64 toolset on the command line.

> `spawn cl.exe ENOENT`

You are probably not building from a command prompt that has your compiler (cl.exe) available. Make sure that you are using **x64 Native Tools Command Prompt for VS 2022** (or run the proper `vcvars.bat` from whatever shell you happen to be using). If the error still shows up, run `cl.exe` from the prompt to check that it is indeed working. If not, your visual studio install might have been set up incorrectly.

