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
   staticruntime "On"
   cppdialect "C++14"
   warnings "Extra"
   buildoptions { "-pedantic" }
   
   filter "configurations:Debug"
      defines { "DEBUG" }
      symbols "On"
      
   filter "configurations:Release"
      flags { "LinkTimeOptimization" }
      defines { "NDEBUG" }
      optimize "Full"
      runtime("Release")
      
      
    -- GCC specific build options.
   configuration "gmake"
        -- Needed on 64-bit platforms to be able
        -- to link static libraries to shared libraries.
        buildoptions { "-fPIC" }
        -- Enables some additional warnings.
        buildoptions { "-Wall" }
        -- Disable some warnings.
        buildoptions 
        { 
            -- Pragma warnings caused by OpenMP support not being enabled.
            "-Wno-unknown-pragmas"
        }