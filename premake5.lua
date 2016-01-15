solution "BigInt"
   configurations { "Debug", "Release" }
   platforms { "x64" }
   flags {"MultiProcessorCompile"}
   location "build/"
   
   
project "BigInt"
   kind "ConsoleApp"
   language "C++"
   targetdir "bin/%{cfg.platform}_%{cfg.buildcfg}"
   objdir "bin/obj/%{cfg.platform}_%{cfg.buildcfg}/%{prj.name}"
   debugdir "bin/%{cfg.platform}_%{cfg.buildcfg}"
   files { "src/**.h", "test/**.cpp" }
   includedirs { "src/" }

   
   filter "configurations:Debug"
      defines { "DEBUG" }
      flags { "Symbols" }
      
   filter "configurations:Release"
      flags { "LinkTimeOptimization", "StaticRuntime" }
      defines { "NDEBUG" }
      optimize "On"
      --linkoptions {"/NODEFAULTLIB:MSVCRT"}
      runtime("Release")
      
      
    -- GCC specific build options.
   configuration "gmake"
        -- Needed on 64-bit platforms to be able
        -- to link static libraries to shared libraries.
        buildoptions { "-fPIC" }
        -- Enables some additional warnings.
        buildoptions { "-Wall" }
        -- Enables C++11 support.
        buildoptions { "-std=c++11" }
        -- Disable some warnings.
        buildoptions 
        { 
            -- Pragma warnings caused by OpenMP support not being enabled.
            "-Wno-unknown-pragmas", 
            -- Comparison between an unsigned and a signed integer.
            "-Wno-sign-compare", 
            -- Conversion between an unsigned and a signed integer.
            "-Wno-sign-conversion",
            -- Unused variables.
            "-Wno-unused-variable",
            -- Unused values.
            "-Wno-unused-value",
            -- Unused but set variable.
            "-Wno-unused-but-set-variable",
            -- Unused functions.
            "-Wno-unused-function",
            -- Breaking strict aliasing rules.
            "-Wno-strict-aliasing",
            -- Compiler warns that it optimizes code based on the 
            -- assumption that signed integer overflows do not occur.
            "-Wno-strict-overflow"
        }