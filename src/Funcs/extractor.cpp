#include "extractor.h"

void print_yellow(const char *str) {
    printf("\033[1;33m%s\033[0m", str);
}

bool extract_boolean_value(const char *json, const char *key, bool *value) {
    char search_str[128];
    snprintf(search_str, sizeof(search_str), "\"%s\":", key);
    const char *key_start = strstr(json, search_str);

    if (!key_start) {
        return false;
    }

    // Extract the line containing the key
    char key_line[128];
    const char *end_of_line = strchr(key_start, '\n');
    if (!end_of_line) {
        end_of_line = key_start + strlen(key_start); // If no newline found, use the end of the string
    }
    size_t len = end_of_line - key_start;
    strncpy(key_line, key_start, len);
    key_line[len] = '\0';

    // Check the boolean value
    char true_str[100], false_str[100];
    snprintf(true_str, sizeof(true_str), "\"%s\": true,", key);
    snprintf(false_str, sizeof(false_str), "\"%s\": false,", key);

    if (strcasecmp(key_line, true_str) == 0) {
        *value = true;
        return true;
    } else if (strcasecmp(key_line, false_str) == 0) {
        *value = false;
        return true;
    } else {
        snprintf(true_str, sizeof(true_str), "\"%s\": true", key);
        snprintf(false_str, sizeof(false_str), "\"%s\": false", key);

        if (strcasecmp(key_line, true_str) == 0) {
            *value = true;
            return true;
        } else if (strcasecmp(key_line, false_str) == 0) {
            *value = false;
            return true;
        }

        return false;
    }
}

bool extract_integer_value(const char *json, const char *key, int *value) {
    char search_str[128];
    snprintf(search_str, sizeof(search_str), "\"%s\":", key);
    const char *key_start = strstr(json, search_str);

    if (!key_start) {
        return false;
    }

    key_start += strlen(search_str);

    if (sscanf(key_start, "%d", value) != 1) {
        return false;
    }

    return true;
}

bool extract_double_value(const char *json, const char *key, double *value) {
    char search_str[128];
    snprintf(search_str, sizeof(search_str), "\"%s\":", key);
    const char *key_start = strstr(json, search_str);

    if (!key_start) {
        return false;
    }

    key_start += strlen(search_str);

    if (sscanf(key_start, "%lf", value) != 1) {
        return false;
    }

    return true;
}

bool extract_string_value(const char *json, const char *key, std::string &name) {
    char search_str[128];
    snprintf(search_str, sizeof(search_str), "\"%s\": \"", key);
    
    const char *key_start = strstr(json, search_str);
    if (!key_start) {
        return false;
    }
    key_start += strlen(search_str);  // Move past the key and the quote

    const char *end_quote = strchr(key_start, '\"');
    if (!end_quote) {
        return false;
    }

    // Calculate the length of the string to extract, not including the quotes
    size_t length = end_quote - key_start;
    
    // Assign the extracted value to 'name', omitting the quotes
    name.assign(key_start, length);

    return true;
}

char *read_file(const char *filename) {
    FILE *file = fopen(filename, "rb"); // Open in binary mode
    if (!file) {
        return NULL;
    }

    if (fseek(file, 0, SEEK_END) != 0) {
        perror("Failed to seek to end of file");
        fclose(file);
        return NULL;
    }

    long length = ftell(file);
    if (length == -1) {
        perror("Failed to tell the position");
        fclose(file);
        return NULL;
    }
    fseek(file, 0, SEEK_SET);

    char *buffer = (char *)malloc(length + 1);
    if (!buffer) {
        perror("Failed to allocate memory");
        fclose(file);
        return NULL;
    }

    size_t result = fread(buffer, 1, length, file);
    if (result != (size_t)length) {
        free(buffer);
        if (feof(file)) {
            perror("Unexpected end of file");
        } else if (ferror(file)) {
            perror("Error reading from file");
        }
        fclose(file);
        return NULL;
    }

    buffer[length] = '\0'; // Null terminate the string

    fclose(file);
    return buffer;
}

char *get_group(const char *json, const char *name) {
    char buffer[512];
    sprintf(buffer, "\"%s\":", name);
    const char *field_ptr = strstr(json, buffer);

    if (!field_ptr) {
        sprintf(buffer, "[Warning]: \"%s\" not found in JSON\n", name);
        print_yellow(buffer);
        return NULL;
    }

    // Find the start of the JSON object
    const char *obj_start = strchr(field_ptr, '{');
    if (!obj_start) {
        print_yellow("[Warning]: Start of JSON object not found\n");
        return NULL;
    }

    // Find the end of the JSON object, which is the next '}' character
    const char *obj_end = strchr(obj_start, '}');
    if (!obj_end) {
        print_yellow("[Warning]: End of JSON object not found\n");
        return NULL;
    }

    // Calculate the length of the JSON object and allocate memory for it
    int obj_length = obj_end - obj_start + 1;
    char *json_object = (char *)malloc(obj_length + 1); // Plus one for the null terminator
    if (!json_object) {
        print_yellow("[Error]: Memory allocation failed\n");
        return NULL;
    }

    // Copy the JSON object to the new buffer
    strncpy(json_object, obj_start, obj_length);
    json_object[obj_length] = '\0'; // Null-terminate the string

    return json_object;
}

int extract_magnetic_field(const char *json) {
    int returnv = 0;

    char *json_object = get_group(json, "magnetic_field");
    if (json_object == NULL) {
        return -1;
    }

    // Use the json_object to extract the values using the extract_double_value function
    if (!extract_double_value(json_object, "Bt", &Bt)) {
        print_yellow("[Warning]: Parameter \"Bt\" not found in \"magnetic field\"\n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Bv", &Bv)) {
        print_yellow("[Warning]: Parameter \"Bv\" not found in \"magnetic field\"\n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Bh", &Bh)) {
        print_yellow("[Warning]: Parameter \"Bh\" not found in \"magnetic field\"\n");
        returnv = -1;
    }

    // Free the allocated memory
    free(json_object);

    return returnv; // Return the result of the extraction process
}

int extract_toroidal_machine_geometry(const char *json) {
    int returnv = 0;

    char *json_object = get_group(json, "toroidal_machine_geometry");

    if (!extract_double_value(json_object, "R", &R)) {
        print_yellow("[Warning]: Parameter \"R\" not found on \"toroidal machine geometry\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "a", &a)) {
        print_yellow("[Warning]: Parameter \"a\" not found on \"toroidal machine geometry\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "b", &b)) {
        print_yellow("[Warning]: Parameter \"b\" not found on \"toroidal machine geometry\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "lHFS", &lHFS)) {
        print_yellow("[Warning]: Parameter \"lHFS\" not found on \"toroidal machine geometry\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "lLFS", &lLFS)) {
        print_yellow("[Warning]: Parameter \"lLFS\" not found on \"toroidal machine geometry\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nlimiters", &nlimiters)) {
        print_yellow("[Warning]: Parameter \"nlimiters\" not found on \"toroidal machine geometry\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Vpl", &Vpl)) {
        print_yellow("[Warning]: Parameter \"Vpl\" not found on \"toroidal machine geometry\" \n");
        returnv = -1;
    }

    free(json_object);

    return returnv; // Success
}

int extract_neutral_pressure(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "neutral_pressure");

    if (!extract_double_value(json_object, "pHe", &pHe)) {
        print_yellow("[Warning]: Parameter \"pHe\" not found on \"neutral pressure\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "pH2", &pH2)) {
        print_yellow("[Warning]: Parameter \"pH2\" not found on \"neutral pressure\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}


int extract_rf_power(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "rf_power");

    if (!extract_double_value(json_object, "Prf", &Prf)) {
        print_yellow("[Warning]: Parameter \"Ptot\" not found on \"rf power\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "freq", &freq)) {
        print_yellow("[Warning]: Parameter \"freq\" not found on \"rf power\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "dtpramp", &dtpramp)) {
        print_yellow("[Warning]: Parameter \"dtpramp\" not found on \"rf power\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bselfcol", &bselfcol)) {
        print_yellow("[Warning]: Parameter \"bselfcol\" not found on \"rf power\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_type(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "type");

    if (!extract_boolean_value(json_object, "bfixpowerfrac", &bfixpowerfrac)) {
        printf("[Warning]: Parameter \"bfixpowerfrac\" not found on \"type\" \n");
        returnv = -1;
    }

    if (!extract_boolean_value(json_object, "bnefix", &bnefix)) {
        printf("[Warning]: Parameter \"bnefix\" not found on \"type\" \n");
        returnv = -1;
    }

    if (!extract_boolean_value(json_object, "bTOMAS", &bTOMAS)) {
        printf("[Warning]: Parameter \"bTOMAS\" not found on \"type\" \n");
        returnv = -1;
    }

    if (!extract_boolean_value(json_object, "bnopower", &bnopower)) {
        printf("[Warning]: Parameter \"bnopower\" not found on \"type\" \n");
        returnv = -1;
    }

    if (!extract_boolean_value(json_object, "bproptone", &bproptone)) {
        printf("[Warning]: Parameter \"bproptone\" not found on \"type\" \n");
        returnv = -1;
    }

    if (!extract_boolean_value(json_object, "bkipt", &bkipt)) {
        printf("[Warning]: Parameter \"bkipt\" not found on \"type\" \n");
        returnv = -1;
    }

    if (!extract_boolean_value(json_object, "blhr", &blhr)) {
        printf("[Warning]: Parameter \"blhr\" not found on \"type\" \n");
        returnv = -1;
    }

    if (!extract_boolean_value(json_object, "bmanuel", &bmanuel)) {
        printf("[Warning]: Parameter \"bmanuel\" not found on \"type\" \n");
        returnv = -1;
    }

    if (!extract_boolean_value(json_object, "bICWC", &bICWC)) {
        printf("[Warning]: Parameter \"bICWC\" not found on \"type\" \n");
        returnv = -1;
    }


    free(json_object);
    return returnv; // Success
}

int extract_general_ec(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "general_ec");

    if (!extract_double_value(json_object, "Rdep", &Rdep)) {
        print_yellow("[Warning]: Parameter \"Rdep\" not found on \"ech\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "pecabs0", &pecabs0)) {
        print_yellow("[Warning]: Parameter \"pecabs0\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "widthech", &widthech)) {
        print_yellow("[Warning]: Parameter \"widthech\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "echbackground", &echbackground)) {
        print_yellow("[Warning]: Parameter \"echbackground\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "harmonic", &harmonic)) {
        print_yellow("[Warning]: Parameter \"harmonic\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "muw", &muw)) {
        print_yellow("[Warning]: Parameter \"muw\" not found on \"ech\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_necfix(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "necfix");

    if (!extract_integer_value(json_object, "ic", &ic)) {
        print_yellow("[Warning]: Parameter \"ic\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "necfix", &necfix)) {
        print_yellow("[Warning]: Parameter \"necfix\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Pini", &Pini)) {
        print_yellow("[Warning]: Parameter \"Pini\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "tauP", &tauP)) {
        print_yellow("[Warning]: Parameter \"tauP\" not found on \"ech\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_tomas(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "tomas");

    if (!extract_double_value(json_object, "Rdep1", &Rdep1)) {
        print_yellow("[Warning]: Parameter \"Rdep1\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "PRdep1", &PRdep1)) {
        print_yellow("[Warning]: Parameter \"PRdep1\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Rdep2", &Rdep2)) {
        print_yellow("[Warning]: Parameter \"Rdep2\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "PRdep2", &PRdep2)) {
        print_yellow("[Warning]: Parameter \"PRdep2\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Rdep3", &Rdep3)) {
        print_yellow("[Warning]: Parameter \"Rdep3\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "PRdep3", &PRdep3)) {
        print_yellow("[Warning]: Parameter \"PRdep3\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Rdep4", &Rdep4)) {
        print_yellow("[Warning]: Parameter \"Rdep4\" not found on \"ech\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "PRdep4", &PRdep4)) {
        print_yellow("[Warning]: Parameter \"PRdep4\" not found on \"ech\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_general_ic(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "general_ic");

    if (!extract_double_value(json_object, "Rant", &Rant)) {
        print_yellow("[Warning]: Parameter \"Rant\" not found on \"ich\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "HtoHD", &HtoHD)) {
        print_yellow("[Warning]: Parameter \"HtoHD\" not found on \"ich\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_at_lhr(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "lhr");

    if (!extract_double_value(json_object, "fracpne", &fracpne)) {
        print_yellow("[Warning]: Parameter \"fracpne\" not found on \"blhr\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "fraclhr", &fraclhr)) {
        print_yellow("[Warning]: Parameter \"fraclhr\" not found on \"blhr\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "widthlhr", &widthlhr)) {
        print_yellow("[Warning]: Parameter \"widthlhr\" not found on \"blhr\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "lhrbackground", &lhrbackground)) {
        print_yellow("[Warning]: Parameter \"lhrbackground\" not found on \"blhr\" \n"); // Key not found
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_other(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "other");

    if (!extract_boolean_value(json_object, "fixedTe", &fixedTe)) {
        print_yellow("[Warning]: Parameter \"fixedTe\" not found on \"other\" \n"); // Key not found
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_diffusion(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "diffusion");

    if (!extract_boolean_value(json_object, "bDfix", &bDfix)) {
        print_yellow("[Warning]: Parameter \"bDfix\" not found on \"diffusion\" \n"); // Key not found
        returnv = -1;
    }
    if (bDfix)
        if (!extract_double_value(json_object, "Dfix", &Dfix)) {
            print_yellow("[Warning]: Parameter \"Dfix\" not found on \"diffusion\" \n"); // Key not found
            returnv = -1;
        }
    if (!extract_boolean_value(json_object, "bDbohm", &bDbohm)) {
        print_yellow("[Warning]: Parameter \"bDbohm\" not found on \"diffusion\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bDscaling", &bDscaling)) {
        print_yellow("[Warning]: Parameter \"bDscaling\" not found on \"diffusion\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Dfact", &Dfact)) {
        print_yellow("[Warning]: Parameter \"Dfact\" not found on \"transport\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_convection(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "convection");

    if (!extract_boolean_value(json_object, "bVfix", &bVfix)) {
        print_yellow("[Warning]: Parameter \"bVfix\" not found on \"diffusion\" \n"); // Key not found
        returnv = -1;
    }
    if (bVfix)
        if (!extract_double_value(json_object, "Vfix", &Vfix)) {
            print_yellow("[Warning]: Parameter \"Vfix\" not found on \"diffusion\" \n"); // Key not found
            returnv = -1;
        }
    if (!extract_boolean_value(json_object, "bVscaling", &bVscaling)) {
        print_yellow("[Warning]: Parameter \"bVscaling\" not found on \"diffusion\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_integer_value(json_object, "veq", &veq)) {
        print_yellow("[Warning]: Parameter \"veq\" not found on \"transport\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Vfact", &Vfact)) {
        print_yellow("[Warning]: Parameter \"Vfact\" not found on \"transport\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_tune_transport(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "tune_d_and_v");

    if (!extract_boolean_value(json_object, "btunedv", &btunedv)) {
        print_yellow("[Warning]: Parameter \"btunedv\" not found on \"transport\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "btunevleft", &btunevleft)) {
        print_yellow("[Warning]: Parameter \"btunevleft\" not found on \"transport\" \n");
        returnv = -1;
    }
    if (!extract_integer_value(json_object, "il", &il)) {
        print_yellow("[Warning]: Parameter \"il\" not found on \"transport\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nelfix", &nelfix)) {
        print_yellow("[Warning]: Parameter \"nelfix\" not found on \"transport\" \n");
        returnv = -1;
    }
    if (!extract_integer_value(json_object, "ir", &ir)) {
        print_yellow("[Warning]: Parameter \"ir\" not found on \"transport\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nerfix", &nerfix)) {
        print_yellow("[Warning]: Parameter \"nerfix\" not found on \"transport\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Vini", &Vini)) {
        print_yellow("[Warning]: Parameter \"Vini\" not found on \"transport\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "tauV", &tauV)) {
        print_yellow("[Warning]: Parameter \"tauV\" not found on \"transport\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Dini", &Dini)) {
        print_yellow("[Warning]: Parameter \"Dini\" not found on \"transport\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "tauD", &tauD)) {
        print_yellow("[Warning]: Parameter \"tauD\" not found on \"transport\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_physics_to_include(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "physics_to_include");

    if (!extract_boolean_value(json_object, "bH", &bH)) {
        print_yellow("[Warning]: Parameter \"bH\" not found on \"physics to include\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bH2", &bH2)) {
        print_yellow("[Warning]: Parameter \"bH2\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bHe", &bHe)) {
        print_yellow("[Warning]: Parameter \"bHe\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bion", &bion)) {
        print_yellow("[Warning]: Parameter \"bion\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bcx", &bcx)) {
        print_yellow("[Warning]: Parameter \"bcx\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "belas", &belas)) {
        print_yellow("[Warning]: Parameter \"belas\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bcoulomb", &bcoulomb)) {
        print_yellow("[Warning]: Parameter \"bcoulomb\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bimpur", &bimpur)) {
        print_yellow("[Warning]: Parameter \"bimpur\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "btranspions", &btranspions)) {
        print_yellow("[Warning]: Parameter \"btranspions\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "btranspneut", &btranspneut)) {
        print_yellow("[Warning]: Parameter \"btranspneut\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bedge", &bedge)) {
        print_yellow("[Warning]: Parameter \"bedge\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bpol", &bpol)) {
        print_yellow("[Warning]: Parameter \"bpol\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "vdrift", &vdrift)) {
        print_yellow("[Warning]: Parameter \"vdrift\" not found on \"physics to include\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bcoll", &bcoll)) {
        print_yellow("[Warning]: Parameter \"bcoll\" not found on \"physics to include\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_initial_conditions(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "initial_conditions");

    if (!extract_double_value(json_object, "rmaxini", &rmaxini)) {
        print_yellow("[Warning]: Parameter \"rmaxini\" not found on \"initial conditions\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "widthini", &widthini)) {
        print_yellow("[Warning]: Parameter \"widthini\" not found on \"initial conditions\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nebackgroundl", &nebackgroundl)) {
        print_yellow("[Warning]: Parameter \"nebackgroundl\" not found on \"initial conditions\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nebackgroundr", &nebackgroundr)) {
        print_yellow("[Warning]: Parameter \"nebackgroundr\" not found on \"initial conditions\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Ta0", &Ta0)) {
        print_yellow("[Warning]: Parameter \"Ta0\" not found on \"initial conditions\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "Te0", &Te0)) {
        print_yellow("[Warning]: Parameter \"Te0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nevac", &nevac)) {
        print_yellow("[Warning]: Parameter \"nevac\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nH0", &nH0)) {
        print_yellow("[Warning]: Parameter \"nH0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nHi0", &nHi0)) {
        print_yellow("[Warning]: Parameter \"nHi0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nH2i0", &nH2i0)) {
        print_yellow("[Warning]: Parameter \"nH2i0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nH3i0", &nH3i0)) {
        print_yellow("[Warning]: Parameter \"nH3i0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nHeII0", &nHeII0)) {
        print_yellow("[Warning]: Parameter \"nHeII0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nHeIII0", &nHeIII0)) {
        print_yellow("[Warning]: Parameter \"nHeIII0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nCII0", &nCII0)) {
        print_yellow("[Warning]: Parameter \"nCII0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nCIII0", &nCIII0)) {
        print_yellow("[Warning]: Parameter \"nCIII0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nCIV0", &nCIV0)) {
        print_yellow("[Warning]: Parameter \"nCIV0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "nCV0", &nCV0)) {
        print_yellow("[Warning]: Parameter \"nCV0\" not found on \"initial conditions\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_edge_conditions(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "edge_conditions");

    if (!extract_double_value(json_object, "RH", &RH)) {
        print_yellow("[Warning]: Parameter \"RH\" not found on \"edge conditions\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "REH", &REH)) {
        print_yellow("[Warning]: Parameter \"REH\" not found on \"edge conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "gEd", &gEd)) {
        print_yellow("[Warning]: Parameter \"gEd\" not found on \"edge conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "gEv", &gEv)) {
        print_yellow("[Warning]: Parameter \"gEv\" not found on \"edge conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "gEdn", &gEdn)) {
        print_yellow("[Warning]: Parameter \"gEdn\" not found on \"edge conditions\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "gEe", &gEe)) {
        print_yellow("[Warning]: Parameter \"gEe\" not found on \"edge conditions\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_simulation_grid(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "simulation_grid");

    if (!extract_integer_value(json_object, "nmeshp", &nmeshp)) {
        print_yellow("[Warning]: Parameter \"nmeshp\" not found on \"number of grid points\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_input_file(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "input_file");

    if (!extract_boolean_value(json_object, "bfinput", &bfinput)) {
        print_yellow("[Warning]: Parameter \"bfinput\" not found on \"input file\" \n"); // Key not found
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_time_step(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "time_step");

    if (!extract_double_value(json_object, "t0", &t0)) {
        print_yellow("[Warning]: Parameter \"t0\" not found on \"time step\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "tmainend", &tmainend)) {
        print_yellow("[Warning]: Parameter \"tmainend\" not found on \"time step\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "accur", &accur)) {
        print_yellow("[Warning]: Parameter \"accur\" not found on \"time step\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "dtmax", &dtmax)) {
        print_yellow("[Warning]: Parameter \"dtmax\" not found on \"time step\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "dtmin", &dtmin)) {
        print_yellow("[Warning]: Parameter \"dtmin\" not found on \"time step\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "dtinit", &dtinit)) {
        print_yellow("[Warning]: Parameter \"dtinit\" not found on \"time step\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_time_step_for_rf_coupling(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "time_step_for_rf_coupling");

    if (!extract_boolean_value(json_object, "bupdateRFstep", &bupdateRFstep)) {
        print_yellow("[Warning]: Parameter \"bupdateRFstep\" not found on \"time step for rf coupling\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_integer_value(json_object, "updateRF", &updateRF)) {
        print_yellow("[Warning]: Parameter \"updateRF\" not found on \"time step for rf coupling\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "dtRF", &dtRF)) {
        print_yellow("[Warning]: Parameter \"dtRF\" not found on \"time step for rf coupling\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "dtRFvar", &dtRFvar)) {
        print_yellow("[Warning]: Parameter \"dtRFvar\" not found on \"time step for rf coupling\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "dtRFmax", &dtRFmax)) {
        print_yellow("[Warning]: Parameter \"dtRFmax\" not found on \"time step for rf coupling\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "dtRFmin", &dtRFmin)) {
        print_yellow("[Warning]: Parameter \"dtRFmin\" not found on \"time step for rf coupling\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "dtRFconv", &dtRFconv)) {
        print_yellow("[Warning]: Parameter \"dtRFconv\" not found on \"time step for rf coupling\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_advanced_time_step_settings(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "advanced_time_step_settings");

    if (!extract_boolean_value(json_object, "dtsmooth", &dtsmooth)) {
        print_yellow("[Warning]: Parameter \"dtsmooth\" not found on \"advanced time step settings\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_double_value(json_object, "shokparam", &shokparam)) {
        print_yellow("[Warning]: Parameter \"shokparam\" not found on \"advanced time step settings\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "maxtstepincrement", &maxtstepincrement)) {
        print_yellow("[Warning]: Parameter \"maxtstepincrement\" not found on \"advanced time step settings\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "btendvar", &btendvar)) {
        print_yellow("[Warning]: Parameter \"btendvar\" not found on \"advanced time step settings\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "convcrit", &convcrit)) {
        print_yellow("[Warning]: Parameter \"convcrit\" not found on \"advanced time step settings\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "convsavetime", &convsavetime)) {
        print_yellow("[Warning]: Parameter \"convsavetime\" not found on \"advanced time step settings\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "baccurvar", &baccurvar)) {
        print_yellow("[Warning]: Parameter \"baccurvar\" not found on \"advanced time step settings\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "accurcrit", &accurcrit)) {
        print_yellow("[Warning]: Parameter \"accurcrit\" not found on \"advanced time step settings\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "minaccur", &minaccur)) {
        print_yellow("[Warning]: Parameter \"minaccur\" not found on \"advanced time step settings\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_output_parameters(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "output_parameters");

    if (!extract_integer_value(json_object, "Nlog", &Nlog)) {
        print_yellow("[Warning]: Parameter \"Nlog\" not found on \"output parameters\" \n"); // Key not found
        returnv = -1;
    }
    if (!extract_integer_value(json_object, "cc", &cc)) {
        print_yellow("[Warning]: Parameter \"cc\" not found on \"output parameters\" \n");
        returnv = -1;
    }
    if (!extract_integer_value(json_object, "Nloopsave", &Nloopsave)) {
        print_yellow("[Warning]: Parameter \"Nloopsave\" not found on \"output parameters\" \n");
        returnv = -1;
    }
    if (!extract_boolean_value(json_object, "bOutdt", &bOutdt)) {
        print_yellow("[Warning]: Parameter \"bOutdt\" not found on \"output parameters\" \n");
        returnv = -1;
    }
    if (!extract_double_value(json_object, "dtsave", &dtsave)) {
        print_yellow("[Warning]: Parameter \"dtsave\" not found on \"output parameters\" \n");
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_solver_parameters(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "solver_parameters");

    if (!extract_double_value(json_object, "solvertolerance", &solvertolerance)) {
        print_yellow("[Warning]: Parameter \"solvertolerance\" not found on \"solver parameters\" \n"); // Key not found
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_input_file_name(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "file_input");


    if (!extract_string_value(json_object, "input_file_path", sinputfile)) {
        print_yellow("[Warning]: Parameter \"input_file_path\" not found on \"file_input\" \n"); // Key not found
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}

int extract_bmanuel_input(const char *json) {
    int returnv = 0;
    char *json_object = get_group(json, "bmanuel_input");

    if (!extract_string_value(json_object, "bmanuel_input_file", smanuel)) {
        print_yellow("[Warning]: Parameter \"bmanuel_input_file\" not found on \"bmanuel_input\" \n"); // Key not found
        returnv = -1;
    }

    free(json_object);
    return returnv; // Success
}
