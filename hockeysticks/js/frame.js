/*
* Url preview script
 * powered by jQuery (http://www.jquery.com)
*
 * written by Alen Grakalic (http://cssglobe.com)
*
 * for more info visit http://cssglobe.com/post/1695/easiest-tooltip-and-image-preview-using-jquery
*
*/
this.screenshotPreview = function(){    
                /* CONFIG */
                               
                                xOffset = 10;
                                yOffset = 30;
                               
                                // these 2 variable determine popup's distance from the cursor
                                // you might want to adjust to get the right result
                               
                /* END CONFIG */
                $("a.screenshot").hover(function(e){
                                this.t = this.title;
                                this.title = "";    
                                var c = (this.t != "") ? "<br/>" + this.t : "";
                                $("body").append("<p id='screenshot'><iframe width ='200px' height='250px' scrolling='no' src='"+ this.rel +"#page=1&toolbar=0&scrollbar=0&navpanes=0&statusbar=0' ></iframe>"+ c +"</p>");                           
                                                                                                               
                                $("#screenshot")
                                                .css("top",(e.pageY - xOffset) + "px")
                                                .css("left",(e.pageX + yOffset) + "px")
                                                .fadeIn("fast");                                                                                
    },
                function(){
                                this.title = this.t;              
                                $("#screenshot").remove();
    });      
                $("a.screenshot").mousemove(function(e){
                                $("#screenshot")
                                                .css("top",(e.pageY - xOffset) + "px")
                                                .css("left",(e.pageX + yOffset) + "px");
                });                                          
};
 
 
// starting the script on page load
$(document).ready(function(){
                screenshotPreview();
});

