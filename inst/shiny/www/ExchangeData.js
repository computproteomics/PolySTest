

// Add a custom message handler to reset file inputs
Shiny.addCustomMessageHandler("resetFileInputHandler", function(id) {      
  const idFile = `#${id}`;
  const idProgress = `${idFile}_progress`;
  const idBar = `${idProgress} .bar`;

  $(idProgress).css("visibility", "hidden");
  $(idBar).css("width", "0%");
  $(idProgress).addClass("active");
  $(idFile).replaceWith(idFile = $(idFile).clone(true));
});



// Add an event listener to read data when it is sent to the app
$(document).on("shiny:connected", function() {
  window.addEventListener("message", displayMessage, false);

  function displayMessage(evt) { 
    // Parse the incoming message
    console.log(evt.data);
    const inmessage = JSON.parse(evt.data);
    console.log(inmessage); 
    console.log("read message");
    
    // Check if the message is a request for results
    if (inmessage === "Retrieve results") {
      // Send a response indicating that the result request was received
      evt.source.postMessage("PolySTest: result request received", evt.origin);
      
      // After a 5 second delay, send the result data
      setTimeout(() => evt.source.postMessage(result_data, evt.origin), 5000);
      
      // Update the Shiny input value to trigger retrieval of the output
      Shiny.setInputValue("retrieve_output", "Get data");
      setTimeout(function() { Shiny.setInputValue("retrieve_output", "Idle") }, 5000);

    } else {
      // Send a response indicating that data was received
      evt.source.postMessage("PolySTest: data received", evt.origin);
      
      // Update the Shiny input value with the received data
      Shiny.setInputValue("extdata", evt.data);
    }
  }
});




 // check ext window and retrieve results
shinyjs.send_results = function(params)
{
  
  result_data = params.dat;
};
 