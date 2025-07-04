#ifndef rotor_h
#define rotor_h

#include <sys/socket.h>
#include <netinet/in.h>

struct Server{
    int domain;
    int service;
    int protocol;
    unsigned long interface;
    int port;
    int backlog;
    
    struct sockaddr_in address;

    int socket;

    void (*launch)(void);

};

#endif