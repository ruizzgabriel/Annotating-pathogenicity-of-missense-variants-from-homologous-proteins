// Search annotation of the variant in gnomAD and ClinVar, this information will be displayed later on the results template
// because this process is slow.

document.addEventListener("DOMContentLoaded", function() {
    const contentDiv = document.getElementById('annotation-content');
    const fetchDataUrl = contentDiv.getAttribute('data-fetch-url');

    console.log('Fetching data from:', fetchDataUrl); // Debugging statement

    fetch(fetchDataUrl)
    .then(response => response.json())
    .then(data => {
            console.log('Data received:', data); // Debugging statement
            const contentDiv = document.getElementById('annotation-content');
            // Clear the loading message
            contentDiv.innerHTML = '';

            // Check if data is present and populate table or message accordingly
            if (data.data) {
                // Build the table if annotation data exists
                let tableHTML = '<table id="pairs_table" class="display"><thead><tr>';
                const columns = Object.keys(data.data[0]);

                // Add table headers
                columns.forEach(col => tableHTML += `<th>${col}</th>`);
                tableHTML += '</tr></thead><tbody>';

                // Add table rows
                data.data.forEach(row => {
                    tableHTML += '<tr>';
                    columns.forEach(col => tableHTML += `<td>${row[col]}</td>`);
                    tableHTML += '</tr>';
                });

                tableHTML += '</tbody></table>';
                contentDiv.innerHTML = tableHTML;
            } else {
                // Display message if no data is available
                contentDiv.innerHTML = `<p>${data.message}</p>`;
            }
        })
        .catch(error => {
            document.getElementById('annotation-content').innerHTML = `
                <p>Error loading data. Please try again later.</p>`;
            console.error('Error fetching annotation data:', error);
        });
});
