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
   files { "src/**.h", "test/**.c" }
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