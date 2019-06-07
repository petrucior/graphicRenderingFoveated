# projectFovea
This project was created to build MMF fovea to graphic rendering using OpenMP and Cuda libraries.

Note: The CUDA implementation has not yet been finalized!

![scene](https://user-images.githubusercontent.com/3810960/59125304-4dec0c80-8938-11e9-8150-e54df91ee6db.png)

Considering 7 levels and fovea positioned in the center of image (0, 0) and level size, W = (30, 30). This scene was 
defined 300 x 300. 

![sceneFoveated](https://user-images.githubusercontent.com/3810960/59125359-8390f580-8938-11e9-8766-91604b9d1df8.png)

time raytracing normal: 0.427374 segundos

time raytracing mmf: 0.193421 segundos

## Package required
- raytracing ( vide https://github.com/petrucior/raytracing.git )

Note: Keep the directories at the root of the operation system ( OS )

## Usage
Compile files
```
make
```
How to execute your file.
```
./main
```
