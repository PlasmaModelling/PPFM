1 Compilare la libreria alglib : 
      da terminale eseguire     
      g++ -c -std=c++20 alglib/source/*.cpp
2 Creare libreria statica :
      ar rcs libalglib.a alglib/source/*.o

>>alglib precompilata il 10/02/2024 e spostata in alglib/

3 Opzioni tasks.json :

"type": "cppbuild",
            "label": "MyBuild G++",
            "command": "/bin/g++",
            "args": [
                "-g",
                "-std=c++20",
                "${workspaceFolder}/*.cpp",
                "-L${workspaceFolder}/alglib",  //include percorso
                "-lalglib",                     //abilita
                "-o",
                "${fileDirname}/rooster"
            ],
//continua... 
