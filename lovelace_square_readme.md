# Lovelace's Square Coding and Data Guidelines

## Overview

This document provides comprehensive guidelines for writing clean, efficient, and maintainable code and datasets for Lovelace's Square. These practices ensure code and data are reusable, well-documented, and support both human collaboration and AI systems like Ada.

## Core Principles

### FAIR Principles
- **Findable**: Clear naming, documentation, and metadata
- **Accessible**: Open source, no hidden folders or password-protected files
- **Interoperable**: Use standard formats (.csv, .py, .xml, .m) and common practices
- **Reusable**: Well-documented, clearly licensed, understandable code

### General Guidelines
- **Consistency**: Use consistent structure and naming across files
- **Clarity**: Choose descriptive names for functions, variables, and files
- **Organization**: Keep related items together in logical folder structures
- **Documentation**: Include explanations for code functionality and data origins
- **Efficiency**: Write clean code without sacrificing readability
- **Open Formats**: Use standard formats (.txt, .csv, .xlsx, .py, .m)

## Code Documentation Standards

### Function Headers Template

```matlab
function [output1, output2] = FunctionName(input1, input2)
% FUNCTIONNAME Brief description of the function
%
% Authors: Your Name
% Date Created: YYYY-MM-DD
% License: Specify your license here
% Version: Specify the version
% Reviewed by Lovelace's Square team: Yes/No
%
% Detailed function description:
% Explain what your function does, algorithms used, and special features.
%
% Args:
%    input1 (type): Description of the first input parameter
%    input2 (type): Description of the second input parameter
%
% Returns:
%    output1 (type): Description of the first output parameter
%    output2 (type): Description of the second output parameter
%
% Example:
%    [result1, result2] = FunctionName(arg1, arg2)
%
% See also: RELATEDFUNCTION1, RELATEDFUNCTION2

% Your code here
end
```

### README Requirements

Every project folder must include a README.txt file containing:
- What the code does
- How to use it
- Installation/setup instructions
- Usage examples
- Contact information

## Code Style Guidelines

### Formatting Rules

#### Indentation
- Use 4 spaces per indentation level
- Maintain consistent indentation throughout

```matlab
% Good
if isValid
    result = processData(data);
else
    result = [];
end

% Bad
if isValid
result = processData(data);
else
result = [];
end
```

#### Line Length
- Keep lines under 75 characters
- Use `...` for line continuation in MATLAB

```matlab
% Good
summaryText = ['The analysis shows a strong correlation ' ...
               'between variables.'];

% Bad
summaryText = 'The analysis shows a strong correlation between variables.';
```

#### Spacing
- Add spaces around operators
- Use blank lines to separate logical sections
- Align related assignments

```matlab
% Good spacing
result = (a + b) / 2;
count  = count + 1;

% Good sectioning
data = load('datafile.mat');

% Normalize the data
normData = (data - mean(data)) / std(data);

% Plot the result
plot(normData);
```

### Naming Conventions

#### Functions and Classes: UpperCamelCase
```matlab
MyAwesomeFunction()
ChemometricAnalyzer()
```

#### Variables and Properties: lowerCamelCase
```matlab
myVariable = 5;
dataMatrix = rand(10);
```

#### Constants: UPPERCASE
```matlab
MAX_ITERATIONS = 100;
DEFAULT_TIMEOUT = 30;
```

#### Descriptive Names
```matlab
% Good
velocity = 12.5;
result = calculateMean(data);

% Bad
x = 12.5;
r = calcM(data);
```

### Code Structure Best Practices

#### Functions vs Scripts
- Use functions instead of scripts for better encapsulation
- One function per file for easier organization
- Use sections (`%%`) in long scripts

```matlab
% Function example (preferred)
function result = CalculateMean(data)
    result = sum(data) / numel(data);
end

% Script with sections
%% Load Data
data = load('datafile.mat');

%% Process Data
normData = (data - mean(data)) / std(data);

%% Plot Results
plot(normData);
```

#### Organization Structure
```
project/
├── README.txt
├── mainScript.m
├── functions/
│   ├── calculateSomething.m
│   └── plotResults.m
└── +preprocessing/
    ├── normalizeData.m
    └── scaleData.m
```

## Version Control Guidelines

- Commit regularly with descriptive messages
- Use branches for features, fixes, or experiments
- Test code before pushing to shared repositories
- Keep main branch clean and stable

## Dataset Guidelines

### Organization Structure
```
data/
├── raw/
│   ├── sample1.csv
│   └── sample2.csv
├── processed/
│   ├── sample1_clean.csv
│   └── sample2_clean.csv
└── docs/
    ├── metadata.json
    └── README.txt
```

### Documentation Requirements

#### README.txt Template
```
# Dataset Title

## Description
Brief description of what the data contains and its purpose.

## Source
- Instrument: [Instrument name]
- Location: [Lab/Institution]
- Date Collected: [YYYY-MM-DD]
- Experiment: [Brief experiment description]

## Authors
- [Author names and affiliations]

## Structure
[Describe folder structure and file contents]

## Format
[Describe file formats and column meanings]

## Usage
[Code examples for loading/using the data]

## License
[License information]

## Contact
[Contact information for questions]
```

#### Metadata.json Template
```json
{
  "filename.csv": {
    "date_collected": "YYYY-MM-DD",
    "instrument": "Instrument name",
    "location": "Lab location",
    "experiment": "Experiment description",
    "columns": {
      "column1": "units or description",
      "column2": "units or description"
    },
    "notes": "Additional notes about the data"
  }
}
```

### Data Preparation
- Use clean, consistent formats (CSV preferred)
- Include clear headers with units
- Remove personal identifiers and sensitive information
- Provide rich metadata explaining collection methods
- Test download links for public accessibility

## Contribution Requirements

### Code Contributions

#### Acceptable Content
- Algorithms and models (dimensionality reduction, regression, classification, clustering)
- Preprocessing techniques (baseline correction, smoothing, normalization)
- Visualization tools (spectra plots, score plots, heatmaps)
- Exploratory tools for pattern discovery
- Validation and evaluation tools
- Helper functions and utilities
- Teaching and learning resources

#### Code Package Structure
```
package_name.zip
├── main_function.m          # Primary code file
├── README.txt              # Detailed documentation
├── example_data.mat        # Sample data (if needed)
└── example_usage.m         # Usage demonstration
```

#### Submission Descriptions

**Short Description**: One concise line summarizing functionality
```
EMSC algorithm for correcting multiplicative/additive effects in spectral data.
```

**Extended Description**: Detailed explanation including:
- Purpose and problem solved
- Algorithm or method used
- Main logic/approach
- Basic usage instructions

### Dataset Contributions

#### Requirements
- Public download link from stable hosting (Zenodo, Figshare, institutional repositories)
- Clean format with logical organization
- Rich metadata and documentation
- Privacy protection (no personal identifiers)
- Usage suggestions and examples

#### Hosting Services
- **Zenodo** (recommended - provides DOI)
- **Figshare**
- **Institutional repositories**
- **Google Drive/Dropbox** (ensure public access)

## Licensing Requirements

### Acceptable Licenses
- All contributions must use recognized open-source or open-data licenses
- Most restrictive acceptable: non-commercial use only (e.g., CC BY-NC 4.0)
- Contributors must have rights to share all submitted content
- Ensure license compatibility for collaborative work

### Common Open Licenses
- **Code**: MIT, Apache 2.0, GPL v3
- **Data**: CC BY 4.0, CC BY-SA 4.0, CC BY-NC 4.0

## Review Process

### Submission Steps
1. Visit submission page on Lovelace's Square
2. Choose "List your code" or "List your dataset"
3. Fill contributor details (authors, affiliation, contact)
4. Provide title and descriptions (short and extended)
5. Specify technical details (license, keywords, version, category)
6. Upload ZIP file (code) or provide download URL (datasets)
7. Submit for review

### Review Criteria
- Completeness and accuracy of documentation
- Code functionality and quality
- Appropriate licensing
- Clear descriptions and examples
- Absence of prohibited content
- File integrity and accessibility

## Quality Assurance

### Code Quality Checklist
- [ ] Functions have proper headers with all required fields
- [ ] Code follows naming conventions
- [ ] Indentation and spacing are consistent
- [ ] Comments explain complex logic
- [ ] Examples demonstrate usage
- [ ] Dependencies are clearly listed
- [ ] Code is organized in logical structure

### Dataset Quality Checklist
- [ ] README.txt file present and complete
- [ ] Metadata file includes all relevant information
- [ ] Data format is clean and consistent
- [ ] Column headers are descriptive with units
- [ ] No sensitive or personal information included
- [ ] Download link tested and publicly accessible
- [ ] License clearly specified

## Best Practices for LLM Integration

### Code Structure for AI Consumption
- Use clear, descriptive function names that indicate purpose
- Include comprehensive docstrings with parameter types
- Provide concrete usage examples
- Structure code in modular, reusable functions
- Use consistent patterns across similar functions

### Documentation for AI Understanding
- Write descriptions in clear, technical language
- Include specific parameter requirements and constraints
- Provide expected input/output formats
- Specify dependencies and requirements explicitly
- Use standard terminology and conventions

### Error Handling and Validation
- Include input validation in functions
- Provide meaningful error messages
- Document expected behavior for edge cases
- Include fallback options where appropriate

## Contact and Support

- **Email**: contact@lovelacesquare.org
- **Questions**: Don't hesitate to ask about licensing, data preparation, or submission concerns
- **Community**: Browse existing contributions for examples of good documentation and formatting

## Impact and Recognition

Contributing to Lovelace's Square:
- Helps the chemometrics community save time and effort
- Advances open science and accessible methods
- Supports education for students and newcomers
- Builds contributor reputation in open science

Remember: Even small contributions can have significant impact. Simple scripts or datasets can save others hours of work or help students understand concepts.