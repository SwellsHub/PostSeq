<!--
*** Thanks for checking out this README Template. If you have a suggestion that would
*** make this better, please fork the repo and create a pull request or simply open
*** an issue with the tag "enhancement".
*** Thanks again! Now go create something AMAZING! :D
-->





<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

[![MIT License][license-shield]][license-url]




<!-- PROJECT LOGO -->
<br />
<p align="center">



  <h3 align="center">PostSeq</h3>

  <p align="center">
    Intuitive Post-Processing for RNAseq!
    <br />
  </p>
</p>



<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Prerequisites](#manual-installation-prerequisites)
  * [Installation](#manual-installation)
* [Usage](#usage)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)



<!-- ABOUT THE PROJECT -->
## About The Project

PostSeq is a webapp hosted at 134.68.250.10:3838/shinyApp which increases the effciency and approachability of RNAseq analysis.
PostSeq allows you to complete differential expression, PCA, pathway analysis and more through a streamlined graphical interface that requires no coding.

Additionally, PostSeq allows you to upload your sequencing data to a MySQL database where it can be retrieved or further analyzed at any time.

### Built With
PostSeq is built in Shiny, a popular package used to build webapps with R code. PostSeq is also linked to a MySQL database to host and retrieve sequencing data
* [RShiny](https://shiny.rstudio.com/)
* [MySQL](https://www.mysql.com/)




<!-- GETTING STARTED -->
## Getting Started

PostSeq can be accessed with no installation at 134.68.250.10:3838/shinyApp . This is the recommended way to use the app and allows you to access all functionalities of the program. 

If you want to host a local instance of PostSeq, the setup process consists of launching a shiny server on your own linux server and populating it with the files in this Repo. Detailed instructions on how to install and maintain a shiny server can be found at [https://docs.rstudio.com/shiny-server/](https://docs.rstudio.com/shiny-server/)

Further instructions for the process of hosting a local instance of PostSeq can be found [here](#manual-installation-prerequisites).


<!-- USAGE EXAMPLES -->
## Usage

PostSeq's main workflow consists of uploading your sequencing data (comma delimited counts file or BAM file), generating the experimental design from this data, choosing the desired forms of analysis, and submitting.

Once you have uploaded a data set to MySQL additional forms of analysis are opened up in the MySQL and Essentiality tabs of the app.

If you are uploading a personal counts file, it needs to be a csv with the first column consisting of gene ids and the first row consisting of sample names. Each subsequent row should contain the raw counts for the correct gene and sample.

If you are uploading BAM files, choose all bam files to be processed after selecting the uploading BAM files checkbox and the 
"choose all BAM files to upload" file input. Please note that due to the size of BAM files, this upload process will take some time. Once Bam files have been uploaded, the generated counts file can be downloaded so the BAM files will not have to be uploaded in future analysis.

After uploading your sequence data, it is time to generate your experimental design (design matrix).
Selecting "Let me autofill sample data" will automatically populate the sample names boxes on the right side of the app with the column names from your counts file. This is recommended as it avoids any naming inconsistencies between the counts and experimental design file. If you wish to input samples manually, modify the "Number of Samples" numeric input and then enter each sample name in the sample names boxes. Note that the sample names in these boxes and your counts file must match **exactly**!




### Manual Installation Prerequisites

This is an example of how to list things you need to use the software and how to install them.
* npm
```sh
npm install npm@latest -g
```

### Manual Installation

1. Get a free API Key at [https://example.com](https://example.com)
2. Clone the repo
```sh
git clone https:://github.com/your_username_/Project-Name.git
```
3. Install NPM packages
```sh
npm install
```
4. Enter your API in `config.js`
```JS
const API_KEY = 'ENTER YOUR API';
```


<!-- CONTRIBUTING -->
## Contributing

Public contributions are the backbone of the open source community. If you want to add something, that's **awesome**.
Here are the instructions on how to do that:

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Spencer Wells - spencerwells@gmail.com

Project Link: [https://github.com/SwellsHub/PostSeq](https://github.com/SwellsHub/PostSeq)



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=flat-square
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
