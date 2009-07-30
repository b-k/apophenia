
/* An experimental means of providing more script-like means of sending
  arguments to a text file.

  See documentation for usage. 

    Copyright (c) 2009 Ben Klemens. Licensed under the modified GNU GPL v2; see COPYING and COPYING2.  
*/

#define apop_varad_head(type, name)  type variadic_##name(variadic_type_##name x)

#define apop_varad_declare(type, name, ...) \
        typedef struct {            \
                    __VA_ARGS__       ;  \
                } variadic_type_##name;     \
    apop_varad_head(type, name);

#define apop_varad_var(name, value) name = x.name ? x.name : (value);
#define apop_varad_link(name,...) variadic_##name((variadic_type_##name) {__VA_ARGS__})


//And for model settings groups:

#define Apop_model_add_group(model, type, ...)  \
    apop_settings_group_alloc(model, #type, type ## _settings_free, type ## _settings_copy, type ##_settings_init ((type ## _settings) {__VA_ARGS__})); 

#define apop_varad_setting(in, out, name, value) out->name = in.name ? in.name : (value);
