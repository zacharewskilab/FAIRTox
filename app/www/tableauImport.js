shinyjs.init = function(){
                var containerDiv = document.getElementById('vizContainer');
                url = 'http://public.tableau.com/views/RegionalSampleWorkbook/Storms';
                var viz = new tableau.Viz(containerDiv, url); }
        