### QuCF

<!--
**QuCF/QuCF** is a ✨ _special_ ✨ repository because its `README.md` (this file) appears on your GitHub profile.

Here are some ideas to get you started:

- 🔭 I’m currently working on ...
- 🌱 I’m currently learning ...
- 👯 I’m looking to collaborate on ...
- 🤔 I’m looking for help with ...
- 💬 Ask me about ...
- 📫 How to reach me: ...
- 😄 Pronouns: ...
- ⚡ Fun fact: ...
-->


QuCF uses a slightly modified version of the [QuEST toolkit](https://github.com/QuEST-Kit/QuEST).

To compile QuCF, create a folder [QuCF]:
- cd [QuCF] 
- git clone git@github.com:QuCF/QuCF.git
- cd QuCF/build_qucf
- source tobuild
- make

REMARK: you might need to correct the `-DGPU_COMPUTE_CAPABILITY` in the `tobuild` files according to your GPU device
(https://developer.nvidia.com/cuda-gpus#compute).

For post-processing QuCF simulations, one also needs [https://github.com/QuCF/scripts-py]:
    - cd [QuCF]
    - git clone https://github.com/QuCF/scripts-py/tree/main
Hence, [QuCF]/scripts-py will contain .py files for post-processing simulations from [QuCF]/QuCF.



