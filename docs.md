# Tiny-Pixels Engine

## Classes and description

### Adjustment (`Adjustment.h`)
Used to store the applicable adjustments to a certain pixel/image.

#### Member variables
The expression after the colon indicates the recommended range. Note that there are limited support for values lying outside the range and should be used with caution. 

- `brightness`: $[-150, 150]$
- `contrast`: $[-150, 150]$
- `hue`: $[-360, 360]$
- `saturation`: $[-1, 1]$
- `value`: $[-1, 1]$
- `lift`: $[-1, 1]$
- `gamma`: $(0, 10]$
- `gain`: $[0.5, 1.5]$

`gamma` and `gain` have default value 1 and all other members default to 0.

#### Member methods
- Constructors

    - `Adjustment();` (Default constructor)

        Initializes an Adjustment object with default values as mentioned above.

    - `Adjustment(double brightness, double contrast, double hue, double saturation, double value, double lift, double gamma, double gain);`

        Initializes an `Adjustment` object with the specified values
- `static Adjustment& create_adj_hsv(double hue, double saturation, double value);`

    Returns an initialized `Adjustment` object with the indicated hue, saturation, and value; and default values for all other data members.

- `static Adjustment& create_adj_bc(double brightness, double contrast);`
    
    Returns an initialized `Adjustment` object with the indicated brightness and contrast values; and default values for all other data members. 

- `static Adjustment& create_adj_lgg(double lift, double gamma, double gain);`
    
    Returns an initialized `Adjustment` object with the indicated lift, gamma, and gain values; and default values for all other data members. 

### Color (`Color.h`)


### Font (`Font.h`)


### Image (`Image.h`)


### Interpolation (`Interpolation.h`)

