#include "GitSHA1.h"

const char* CMSat::get_version_sha1()
{
    static const char myversion_sha1[] = "81464784c3ffc64d510488ed9d2ba85935bd0312";
    return myversion_sha1;
}

const char* CMSat::get_version_tag()
{
    static const char myversion_tag[] = "5.11.11";
    return myversion_tag;
}

const char* CMSat::get_compilation_env()
{
    return " ";
}

