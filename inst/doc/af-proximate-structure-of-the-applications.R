## -----------------------------------------------------------------------------
#| label: fig-filestructure
#| out-width: "80%"
#| fig-cap: "Files inside the application file"
#| echo: false
#| fig-align: center
#| fig-retina: 0.85
#| out-extra: 'style="background-color: #FFFFFF; border: 10px solid transparent; padding:0px; display: inline-block;"'

knitr::include_graphics("file_structure.jpg")


## -----------------------------------------------------------------------------
#| eval: true

first_pixel_nir <- 4
first_pixel_count_nir <- first_pixel_nir + 1
first_wavelength_nir <- 2.04E-10 * first_pixel_count_nir^5 + 
  -1.28E-07 * first_pixel_count_nir^4 + 
   2.80E-05 * first_pixel_count_nir^3 +
  -4.76E-03 * first_pixel_count_nir^2 + 
   3.89     * first_pixel_count_nir   + 
   880.06
first_wavelength_nir

first_pixel_count_vis <- 823
first_wavelength_vis <- 0.0 * first_pixel_count_vis^4 + 
   0.0          * first_pixel_count_vis^3 +
  -7.586146E-05 * first_pixel_count_vis^2 + 
   2.12726      * first_pixel_count_vis   + 
  -1301.079
first_wavelength_vis


## -----------------------------------------------------------------------------
#| eval: true

pixel_sequence_nir <- 4:272
pixel_sequence_count_nir <- pixel_sequence_nir + 1

wavelengths_nir <- 2.04E-10 * pixel_sequence_count_nir^5 + 
  -1.28E-07 * pixel_sequence_count_nir^4 +
   2.80E-05 * pixel_sequence_count_nir^3 +
  -4.76E-03 * pixel_sequence_count_nir^2 + 
   3.89     * pixel_sequence_count_nir   + 
  880.06

wavelengths_nir[1:5]
wavelengths_nir[(length(wavelengths_nir) - 5):length(wavelengths_nir)]

pixel_sequence_vis <- 823:1074
wavelengths_vis <- 0.0 * pixel_sequence_vis^4 + 
   0.0          * pixel_sequence_vis^3 +
  -7.586146E-05 * pixel_sequence_vis^2 + 
   2.12726      * pixel_sequence_vis   + 
  -1301.079

wavelengths_vis[1:5]
wavelengths_vis[(length(wavelengths_vis) - 5):length(wavelengths_vis)]

