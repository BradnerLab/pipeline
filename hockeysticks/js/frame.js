



this.screenshotPreview = function(filename){    
                               
    xOffset = 10;
    yOffset = 30;
       // these 2 variable determine popup's distance from the cursor
    console.log("here")


    $("#linegraph").append("<p><embed src=" + filename + " type='application/pdf'></p>");                           
                                                                                   
    // $("#screenshot")
    //     .css("top",(e.pageY - xOffset) + "px")
    //     .css("left",(e.pageX + yOffset) + "px")
    //     .fadeIn("fast");                                                                                                                       
};
 
 
// // starting the script on page load
// $(document).ready(function(){
//     screenshotPreview();
// });

