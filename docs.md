# Tiny-Pixels Engine

## Classes and description

### Adjustment (`Adjustment.h`)
Stores the applicable adjustments to a certain pixel/image.

#### Member variables
The expression after the colon indicates the recommended range. Note that there are limited support for values lying outside the range and should be used with caution. 

- `double brightness`: $[-150, 150]$
- `double contrast`: $[-150, 150]$
- `double hue`: $[-360, 360]$
- `double saturation`: $[-1, 1]$
- `double value`: $[-1, 1]$
- `double lift`: $[-1, 1]$
- `double gamma`: $(0, 10]$
- `double gain`: $[0.5, 1.5]$

`gamma` and `gain` have default value 1 and all other members default to 0.

#### Member methods
- Constructors

    - `Adjustment();` (Default constructor)

        Initializes an `Adjustment` object with default values as mentioned above.

    - `Adjustment(double brightness, double contrast, double hue, double saturation, double value, double lift, double gamma, double gain);`

        Initializes an `Adjustment` object with the specified values
- `static Adjustment& create_adj_hsv(double hue, double saturation, double value);`

    Returns an initialized `Adjustment` object with the indicated hue, saturation, and value; and default values for all other data members.

- `static Adjustment& create_adj_bc(double brightness, double contrast);`
    
    Returns an initialized `Adjustment` object with the indicated brightness and contrast values; and default values for all other data members. 

- `static Adjustment& create_adj_lgg(double lift, double gamma, double gain);`
    
    Returns an initialized `Adjustment` object with the indicated lift, gamma, and gain values; and default values for all other data members. 

### Color (`Color.h`)
Stores data for a particular color in RGBA format. 

#### Member variables

Representing the RGBA channels of a certain color:

- `double r`
- `double g`
- `double b`
- `double a`

Note that it may also represent color in HSV where `r`, `g`, `b` corresponds to hue, saturation, and value, with recommended range $[0, 360)$, $[0, 1]$, and $[0, 1]$ respectively. 

#### Member methods

- Constructors

    - `Color();` (Default constructor)

        Initializes a `Color` object with default values of `r = 0`, `g = 0`, `b = 0`, and `a = 255`.

    - `Color(double r, double g, double b);`

        Initializes a `Color` object with the specified RGB values, `a` defaults to 255.

    - `Color(double r, double g, double b, double a);` 

        Initializes a `Color` object with the specified RGBA values. 
        
    - `Color(Color* color)`

        Initializes a copy of the provided `color`. 

- `void hsv_to_rgb(double h_deg, double s, double v);`
    
    Alters the calling object to be the RGB equivalent of the provided HSV values. 

- `Color& to_hsv();`
    
    Alters the calling object to be the HSV equivalent of the currently stored RGB values. 

- `Color& rgb_to_hsv(double rr, double gg, double bb);`

    Returns a new `Color` object that stores the HSV equivalent of the provided RGB values. 

- `double get(int col);`

    A numeric interface to retrieving member variables with the mapping 

    - 0 $\rightarrow$ `r`
    - 1 $\rightarrow$ `g`
    - 2 $\rightarrow$ `b`
    - 3 $\rightarrow$ `a`

    All other values default to -1

- `void set(double rr, double gg = 0, double bb = 0, double aa = 255);`
    
    Sets the calling object's member variables to the provided values.

- `void set(Color color);`

    Sets the calling object's member variables to match that of the provided `color`.

- `bool set(int col, double val);`

    A numeric interface similar to `Color::get(int col)` that sets the specified variable to the provided value. 

    Return value indicates if the operation was successful.

- `bool operator<(const Color& other) const;`

    Comparison operator that compares the calling object with the provided `color` by `r`, `g`, `b`, in that order. 

- `Color operator+(const Color& other) const;`
    
    Adds two colors component wise.

- `Color operator*(const Color& other) const;`

    Multiplies two colors component wise.

- `Color operator*(double mult) const;`

    Multiplies all member variables in the calling object by `mult`.

- `Color operator/(const Color& other) const;`

    Divides two colors component wise.

- `Color operator/(double div) const;`

    Divides all member variables in the calling object by `div`.

- `static double luminance(double r, double g, double b);`

    Calculates luminance from the given RGB values via $0.2126r + 0.7152g + 0.0722b$.

- `double luminance();`

    Calculates luminance of the calling object.

- `static double luminance(const Color& c);`

    Calculates luminance of the provided `Color` object.

- `Color& apply_adj_rgb(Adjustment adj, double factor);`

    Applies the provided adjustment to the calling color accounting for the provided factor, with 1 being fully applied, and 0 being none applied. 

### Font (`Font.h`)

#### Member variables

Font object from the `schrift` library.

- `SFT sft = {NULL, 12, 12, 0, 0, SFT_DOWNWARD_Y | SFT_RENDER_IMAGE};`

#### Member functions

- Constructor
    
    - `Font(const char* fontfile, uint16_t size);`

        Initializes the font object with the provided TrueType Font (`.ttf`) with the provided size.

- `void setSize(uint16_t size);`

    Sets the font size of the calling object to the provided size.

### Image (`Image.h`)



### Interpolation (`Interpolation.h`)

