# Epigenetic State of a Fatty Acid Synthesis (FAS) Machinery in Glioblastoma Multiforme (GBM)

## Introduction

Altered lipid metabolism is an emerging hallmark of cancer. In GBM, lipid metabolism is essential to cancer cell proliferation. Using two independent population-scale GM cohorts, this study stratifies each into _predominant subgroups R and L_ by applying unsupervised machine learning to the epigenetic profile of a high-profile candidate gene in FAS. 

Subsequently, this study compared the FAS epigenetic subgroups, one cohort at a time, and...

* Shows significant differences in _cancer stemness_
* Shows significant differences in _cellular differentiation potential_ 
* Identifies distinct _neuronal cell type proportions_
* Demonstrates a significant overlap between subgrouping and global epigenetic landscape 
* Uses EWAS to reveal biologically relevant gene sets and processes 
* Reveals differences in _overall survival_ by subgroups

## Repository Structure

The scripts prefixed with numbers, e.g. `01_`, are used for the analyses presented in the final write-up. File `utils.R` is loaded as a utility module with shared functions and global variables per R session. Some variables in this script are masked.

Scripts for data cleaning, wrangling, exploratory data analyses (EDAs), analyses not presented in the final write-up, and additional data/results I/O are not included.

## Audience Engagement

### Preprint and Manuscript

This work has been submitted for publication and under review. A preprint version of this work is available at [Preprints](https://www.preprints.org/manuscript/202509.0713).

### Interactive Dashboard

An interactive dashboard of crucial measures, analogous to "KPIs" in business settings, can be found in this [Tableau Public dashboard](https://public.tableau.com/app/profile/david.chen8785/viz/efas/Final):


<html>
<div class='tableauPlaceholder' id='viz1770393353259' style='position: relative'>
    <noscript>
        <a href='#'><img alt='Final ' src='https:&#47;&#47;public.tableau.com&#47;static&#47;images&#47;ef&#47;efas&#47;Final&#47;1_rss.png' style='border: none' /></a>
    </noscript>
    <object class='tableauViz' style='display:none;'>
        <param name='host_url' value='https%3A%2F%2Fpublic.tableau.com%2F' />
        <param name='embed_code_version' value='3' />
        <param name='site_root' value='' />
        <param name='name' value='efas&#47;Final' />
        <param name='tabs' value='no' />
        <param name='toolbar' value='yes' />
        <param name='static_image' value='https:&#47;&#47;public.tableau.com&#47;static&#47;images&#47;ef&#47;efas&#47;Final&#47;1.png' />
        <param name='animate_transition' value='yes' />
        <param name='display_static_image' value='yes' />
        <param name='display_spinner' value='yes' />
        <param name='display_overlay' value='yes' />
        <param name='display_count' value='yes' />
        <param name='language' value='en-US' />
        <param name='filter' value='publish=yes' />
    </object>
</div>
<script type='text/javascript'>
    var divElement = document.getElementById('viz1770393353259');                    var vizElement = divElement.getElementsByTagName('object')[0];                    if ( divElement.offsetWidth > 800 ) { vizElement.style.width='100%';vizElement.style.height=(divElement.offsetWidth*0.75)+'px';} else if ( divElement.offsetWidth > 500 ) { vizElement.style.width='100%';vizElement.style.height=(divElement.offsetWidth*0.75)+'px';} else { vizElement.style.width='100%';vizElement.style.height='2427px';}                     var scriptElement = document.createElement('script');                    scriptElement.src = 'https://public.tableau.com/javascripts/api/viz_v1.js';                    vizElement.parentNode.insertBefore(scriptElement, vizElement);
</script>
</html>

### Responsible Use

Use the data/results, findings, and interpretations of this work responsibly and at your own risk. Cite relevant sources wherevera applicable.

