



this.screenshotPreview = function(filename){    
                               
    xOffset = 10;
    yOffset = 30;
       // these 2 variable determine popup's distance from the cursor
    console.log("here")

    $("#linegraph").append("<p id='screenshot'><iframe width ='200px' height='250px' scrolling='no' src='"+ filename +"#page=1&toolbar=0&scrollbar=0&navpanes=0&statusbar=0' ></iframe>" +"</p>");                           
                                                                                   
    $("#screenshot")
        .css("top",(e.pageY - xOffset) + "px")
        .css("left",(e.pageX + yOffset) + "px")
        .fadeIn("fast");                                                                                                                       
};
 
 
// // starting the script on page load
// $(document).ready(function(){
//     screenshotPreview();
// });

