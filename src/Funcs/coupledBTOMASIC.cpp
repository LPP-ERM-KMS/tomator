#include "coupledpower.h"
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

void bTOMASIC_func() {
    //# create data pipes for inter-process communication (IPC) (comm via RAM instead of disk):
    //## Assign pipe filenames in tmp 
    char RaPipe[] = "/tmp/RaNamedPipe";
    char NePipe[] = "/tmp/NeNamedPipe";
    char TePipe[] = "/tmp/TeNamedPipe";
    char NiPipe[] = "/tmp/NiNamedPipe";
    char TiPipe[] = "/tmp/TiNamedPipe";
    char PRFePipe[] = "/tmp/PRFeNamedPipe";
    char PRFH2iPipe[] = "/tmp/PRFH2iNamedPipe";
    //## mkfifo(<pathname>, <permission>), (First In, First Out)
    mkfifo(RaPipe, 0777);
    mkfifo(NePipe, 0777);
    mkfifo(TePipe, 0777);
    mkfifo(NiPipe, 0777);
    mkfifo(TiPipe, 0777);
    mkfifo(PRFePipe, 0777);
    mkfifo(PRFH2iPipe, 0777);
    //# Write Density and Temperature to data pipes and after read power depo
    //## Open FIFOs for write only (TW=To Write)
    int RaTW;
    int NeTW;
    int TeTW;
    int NiTW;
    int TiTW;
    RaTW = open(RaPipe, O_RDWR);
    NeTW = open(NePipe, O_RDWR);
    TeTW = open(TePipe, O_RDWR);
    NiTW = open(NiPipe, O_RDWR);
    TiTW = open(TiPipe, O_RDWR);
    //## Open FIFOs for read only (TR=To Read)
    int PRFeTR;
    int PRFH2iTR;
    PRFeTR = open(PRFePipe, O_RDWR);
    PRFH2iTR = open(PRFH2iPipe, O_RDWR);
    //## Convert values to chars
    //### Assign chars
    char RaChar[sizeof(aR)];
    char NeChar[sizeof(nr.ne)];
    char TeChar[sizeof(Tr.Te)];
    char NiChar[sizeof(nr.nH2i)];
    char TiChar[sizeof(Tr.TH2i)];
    memcpy(RaChar,&aR,sizeof(aR));
    memcpy(NeChar,&nr.ne,sizeof(nr.ne));
    memcpy(TeChar,&Tr.Te,sizeof(Tr.Te));
    memcpy(NiChar,&nr.nH2i,sizeof(nr.nH2i));
    memcpy(TiChar,&Tr.TH2i,sizeof(Tr.TH2i));
    //## Write the current values to the FIFOs
    //## write(pipe,buffer,num_bytes)
    write(RaTW, RaChar, strlen(RaChar)+1);
    write(NeTW, NeChar, strlen(NeChar)+1);
    write(TeTW, TeChar, strlen(TeChar)+1);
    write(NiTW, NiChar, strlen(NiChar)+1);
    write(TiTW, TiChar, strlen(TiChar)+1);
    //# Execute ngsolve sim and write result, this script reads the pipes
    //## Execution 
    cout << "Executing python" << endl;
    std::string ScriptPosBase = std::getenv("TOMATORSOURCE");
    std::string ScriptPosADD = "/src/ICSim/TOMASH.py";
    std::string ScriptPos = ScriptPosBase + ScriptPosADD;
    std::string command = "python ";
    command += ScriptPos;
    system(command.c_str());
    //## Close the write pipes
    close(RaTW);
    close(NeTW);
    close(TeTW);
    close(NiTW);
    close(TiTW);
    //## Python script writes to pipes, reading result:
    char PRFe_array_char[80];
    char PRFH2i_array_char[80];
    cout << "Reading Pipes" << endl;
    read(PRFeTR, PRFe_array_char , sizeof(PRFe_array_char));
    read(PRFH2iTR, PRFH2i_array_char , sizeof(PRFH2i_array_char));
    //## Test if it works
    cout << "Printing output" << endl;
    printf(PRFe_array_char);
    close(PRFeTR);
    close(PRFH2iTR);
}
